function [params, robot] = get_inertial_params(robot, varargin)
    % author: jon michaux
    % updated: 20220322 by patrick holmes
    % -- generalizing this function to work with different set types
    % -- specifying uncertainty for individual links
    
    %% parse inputs
    default_set_type = 'point'; % choose 'point', 'interval', or 'polynomial_zonotope'
    default_add_uncertainty_to = 'none'; % choose 'none', 'all', or 'link'
    default_links_with_uncertainty = {}; % if adding uncertainty to specific links, specify names here
    default_mass_range = [1 1];
    default_com_range = [1 1];
    default_track_inertial_generators = false; % PZ specific, either treat inertial gens as dependent or independent
    default_zono_order = 40; % PZ specific, choose number of generators (order * dimension = max # of generators)
    
    p = inputParser;
    addParameter(p, 'set_type', default_set_type, @ischar);
    addParameter(p, 'add_uncertainty_to', default_add_uncertainty_to, @ischar);
    addParameter(p, 'links_with_uncertainty', default_links_with_uncertainty, @(C) all(cellfun(@ischar, C)));
    addParameter(p, 'mass_range', default_mass_range, @(x) all(size(x) == [1, 2]));
    addParameter(p, 'com_range', default_com_range, @(x) all(size(x) == [1, 2]));
    addParameter(p, 'track_inertial_generators', default_track_inertial_generators, @islogical);
    addParameter(p, 'zono_order', default_zono_order, @isscalar);
    parse(p, varargin{:});
    
    % some sanity checks...
    if strcmp(p.Results.set_type, 'point') && ((p.Results.mass_range(1) ~= p.Results.mass_range(2)) || p.Results.com_range(1) ~= p.Results.com_range(2))
        error('If the set type is a point, first and second elements of mass_range and com_range must be identical.')
    end
    if strcmp(p.Results.add_uncertainty_to, 'link') && isempty(p.Results.links_with_uncertainty)
        error('If adding uncertainty to specific links, specify the names of the bodies.');
    end
    if (p.Results.com_range(1) ~= 1 || p.Results.com_range(2) ~= 1)
        warning('Need to double check that this is implemented correctly when moving COM location');
    end
    
    robot.DataFormat = 'col';
    robot.Gravity = [0 0 -9.81];

    num_joints = length(robot.Bodies);
     
    %% initialize placeholders
    
    if strcmp(p.Results.set_type, 'polynomial_zonotope')
        id = [];
        id_names = [];
    end
    
    % moments of inertia
    I = {};
    
    % for spatial
    G = {};
    
    %% get/compute params
    
    for i = 1:num_joints
        %~~~~~~~~~~~~~~~~~inertial parameters
        switch p.Results.set_type
            case 'point'
                if strcmp(p.Results.add_uncertainty_to, 'all') || (strcmp(p.Results.add_uncertainty_to, 'link') && any(strcmp(p.Results.links_with_uncertainty, robot.Bodies{i}.Name)))
                    % mass
                    mass(i) = p.Results.mass_range(1)*robot.Bodies{i}.Mass;

                    % center of mass
                    com(:, i) = p.Results.com_range(1)*robot.Bodies{i}.CenterOfMass';

                    % inertia
                    I_tmp = parallel_axis(robot.Bodies{i}.Inertia, robot.Bodies{i}.Mass, robot.Bodies{i}.CenterOfMass');  
                    I{i} = p.Results.mass_range(1)*I_tmp;  

                    % spatial inertias
                    G{i} = [   I{i}        zeros(3,3); 
                            zeros(3,3)   mass(i)*eye(3)];
                else
                    % mass
                    mass(i) = robot.Bodies{i}.Mass;

                    % center of mass
                    com(:, i) = robot.Bodies{i}.CenterOfMass';

                    % inertia
                    I{i} = parallel_axis(robot.Bodies{i}.Inertia, mass(i), com(:,i));  

                    % spatial inertias
                    G{i} = [   I{i}        zeros(3,3); 
                            zeros(3,3)   mass(i)*eye(3)];
                end
            case 'interval'
                if strcmp(p.Results.add_uncertainty_to, 'all') || (strcmp(p.Results.add_uncertainty_to, 'link') && any(strcmp(p.Results.links_with_uncertainty, robot.Bodies{i}.Name)))
                    [mass_tmp, com_tmp, I_tmp] = get_uncertain_interval_params(robot, p, i);
                    mass(i) = mass_tmp;
                    com(:, i) = com_tmp;
                    I{i} = I_tmp;

                    % spatial inertias
                    G{i} = [   I{i}        zeros(3,3); 
                            zeros(3,3)   mass(i)*eye(3)];
                else
                    % mass
                    mass(i) = interval(robot.Bodies{i}.Mass, robot.Bodies{i}.Mass);

                    % center of mass
                    com(:, i) = interval(robot.Bodies{i}.CenterOfMass', robot.Bodies{i}.CenterOfMass');

                    % inertia
                    I_tmp = parallel_axis(robot.Bodies{i}.Inertia, robot.Bodies{i}.Mass, robot.Bodies{i}.CenterOfMass');
                    I{i} = interval(I_tmp, I_tmp);
                    
                    % spatial inertias
                    G{i} = [   I{i}        zeros(3,3); 
                            zeros(3,3)   mass(i)*eye(3)];
                end
            case 'polynomial_zonotope'
                if strcmp(p.Results.add_uncertainty_to, 'all') || (strcmp(p.Results.add_uncertainty_to, 'link') && any(strcmp(p.Results.links_with_uncertainty, robot.Bodies{i}.Name)))
                   % very similar to interval params, but cast as PZ
                    [mass_tmp, com_tmp, I_tmp] = get_uncertain_interval_params(robot, p, i);
                    dim = 3;
                    
                    % center of mass for each link
                    inf = com_tmp.inf;
                    sup = com_tmp.sup;
                    c = (inf+sup)./2;
                    g = [];
                    expMat = [];
                    for j = 1:dim
                        if sup(j) ~= inf(j)
                            g(:, j) = zeros(dim, 1);
                            g(j, j) = (sup(j) - inf(j))./2;
                        end
                    end
                    
                    if p.Results.track_inertial_generators
                        [id, id_names, id_add] = update_id(id, id_names, size(g, 2), ['c', num2str(i)]);
                        com{i, 1} = polyZonotope_ROAHM(c, g, [], eye(size(g, 2)), id_add);
                    else
                        com{i, 1} = polyZonotope_ROAHM(c, [], g);
                    end
                    
                    % link masses
                    inf = mass_tmp.inf;
                    sup = mass_tmp.sup;
                    C = (inf + sup)./2 * eye(dim);
                    G = [];
                    if sup ~= inf
                        G = cat(3, G, (sup - inf)./2*eye(dim));
                    end
                    if p.Results.track_inertial_generators
                        [id, id_names, id_add] = update_id(id, id_names, size(G, 3), ['m', num2str(i)]);
                        mass{i, 1} = matPolyZonotope_ROAHM(C, G, [], eye(size(G, 3)), id_add);
                    else
                        mass{i, 1} = matPolyZonotope_ROAHM(C, [], G);
                    end
                    
                    % inertia wrt to center of mass frame
                    inf = I_tmp.inf;
                    sup = I_tmp.sup;
                    C = (inf + sup)./2;
                    G = [];
                    for j = 1:dim
                        for k = j:dim
                            if sup(j,k) ~= inf(j, k)
                                G_tmp = zeros(dim, dim);
                                G_tmp(j, k) = (sup(j,k) - inf(j,k))./2;
                                if j ~= k
                                    G_tmp(k, j) = G_tmp(j, k);
                                end
                                G = cat(3, G, G_tmp);
                            end
                        end
                    end
                    if p.Results.track_inertial_generators
                        [id, id_names, id_add] = update_id(id, id_names, size(G, 3), ['I', num2str(i)]);
                        I{i, 1} = matPolyZonotope_ROAHM(C, G, [], eye(size(G, 3)), id_add);
                    else
                        I{i, 1} = matPolyZonotope_ROAHM(C, [], G);
                    end
                else
                    % use nominal parameters as centers
                    mass{i, 1} = matPolyZonotope_ROAHM(robot.Bodies{i}.Mass*eye(3));
                    com{i, 1} = polyZonotope_ROAHM(robot.Bodies{i}.CenterOfMass');
                    I_tmp = parallel_axis(robot.Bodies{i}.Inertia, robot.Bodies{i}.Mass, robot.Bodies{i}.CenterOfMass');
                    I{i, 1} = matPolyZonotope_ROAHM(I_tmp);
                end
            otherwise
                error('Set type not recognized.');
        end
                       
    end
    
    params.set_type = p.Results.set_type;
    params.mass = mass;
    params.com = com;
    params.I = I;
    params.G = G;
    params.add_uncertainty_to = p.Results.add_uncertainty_to;
    params.links_with_uncertainty = p.Results.links_with_uncertainty;
    params.mass_range = p.Results.mass_range;
    params.com_range = p.Results.com_range;
    params.zono_order = p.Results.zono_order;
    if strcmp(p.Results.set_type, 'polynomial_zonotope')
        params.id = id;
        params.id_names = id_names;
    end
end

function [mass, com, I] = get_uncertain_interval_params(robot, p, i)
    % interval mass
    mass = interval(p.Results.mass_range(1), p.Results.mass_range(2))*robot.Bodies{i}.Mass;
    
    % interval COM
    com = interval(p.Results.com_range(1), p.Results.com_range(2))*robot.Bodies{i}.CenterOfMass';
    
    % interval inertia
    % pat 20220323 explanation for changes: 
    % here I'm assuming the mass distribution about the COM is fixed,
    % (except the mass is not exactly known). so the inertia matrix does
    % not depend on the COM position. uncertainty in the COM position
    % corresponds to the link's COM not being exactly known relative to the
    % joint position.
    I_tmp = parallel_axis(robot.Bodies{i}.Inertia, robot.Bodies{i}.Mass, robot.Bodies{i}.CenterOfMass');  
    I = interval(p.Results.mass_range(1), p.Results.mass_range(2))*I_tmp;
    
end