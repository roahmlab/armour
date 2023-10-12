function [params, robot] = get_robot_params_v2(robot, varargin)
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
    default_zono_order = 40; % PZ specific, choose number of generators (order * dimension = max # of generators)
    
    p = inputParser;
    addParameter(p, 'set_type', default_set_type, @ischar);
    addParameter(p, 'add_uncertainty_to', default_add_uncertainty_to, @ischar);
    addParameter(p, 'links_with_uncertainty', default_links_with_uncertainty, @(C) all(cellfun(@ischar, C)));
    addParameter(p, 'mass_range', default_mass_range, @(x) all(size(x) == [1, 2]));
    addParameter(p, 'com_range', default_com_range, @(x) all(size(x) == [1, 2]));
    addParameter(p, 'zono_order', default_zono_order, @isscalar);
    parse(p, varargin{:});
    
    % some sanity checks...
    if strcmp(p.Results.set_type, 'point') && ((p.Results.mass_range(1) ~= p.Results.mass_range(2)) || p.Results.com_range(1) ~= p.Results.com_range(2))
        error('If the set type is a point, first and second elements of mass_range and com_range must be identical.')
    end
    if strcmp(p.Results.set_type, 'point') && ~strcmp(p.Results.add_uncertainty_to, 'none')
        error('Can''t add uncertainty when set_type is point');
    end
    if strcmp(p.Results.add_uncertainty_to, 'link') && isempty(p.Results.links_with_uncertainty)
        error('If adding uncertainty to specific links, specify the names of the bodies.');
    end
    
%     robot = importrobot(robot_name);
    robot.DataFormat = 'col';
    robot.Gravity = [0 0 -9.81];

    %% get joint types
%     [active_joints, fixed_joints] = get_joint_types(robot);
%     num_joints = length(active_joints);
    num_joints = length(robot.Bodies);
    joint_types = cell(num_joints, 1);
    q_index = nan(num_joints, 1);
    q_count = 0;
    for i = 1:robot.NumBodies
       joint_types{i, 1} = robot.Bodies{i}.Joint.Type;
       if ~strcmp(robot.Bodies{i}.Joint.Type, 'fixed')
           % get corresponding configuration for this index.
           q_count = q_count + 1;
           q_index(i, 1) = q_count;
       end
    end
    num_q = q_count;
     
    %% initialize placeholders
    params = {};
%     mass = zeros(1, num_joints);
%     com = zeros(3, num_joints);
    joint_axes = zeros(3, num_joints);
    
    % transforms when robot is in home configuration
    H = {};
    
    % position of the origin of frame {i} with respect to frame {i-1}
    P = zeros(3, num_joints + 1);
    
    % spatial inertias
    G = {};
    
    % moments of inertia
    I= {};
    
    % orientation of com frame {i} in com frame {i-1}
    M = repmat(eye(4), [1,1,num_joints+1]);
    
    % screw axes in world frame
    S = zeros(6, num_joints);
    
    % orientation joint frame{i} in space frame {0}
    Tj = zeros(4, 4, num_joints);
    
    % orientation of com frame {i} in space frame {0}
    Tc = zeros(4,4,num_joints);
    
    %% get/compute params
    
    q0 = homeConfiguration(robot);
    T0 = zeros(4,4,num_joints);
    
    for i = 1:num_joints
        %~~~~~~~~~~~~~~~~~transforms 
        % calculate frame-to-frame transforms
        if i == 1
            H{i} = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Base.Name);      
        else
            H{i} = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Bodies{i-1}.Name);
        end   
        
        T0(:, :, i) = H{i};
        
%         if i == 1
%             T0(:,:,i) = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Base.Name);
%         else
%             T0(:,:,i) = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Bodies{i-1}.Name);
%         end
        
        % position of the origin of frame {i} with respect to frame {i-1}
        P(:, i) = H{i}(1:3, 4);
        
        % axis of rotation
        joint_axes(:,i) = robot.Bodies{i}.Joint.JointAxis';
             
        % orientation joint frame{i} in space frame {0}
        Tj(:,:,i) = getTransform(robot, zeros(num_q,1), robot.Bodies{i}.Name, robot.Base.Name);      
        
        % orientation of com frame {i} in space frame {0}
        K = eye(4);
        K(1:3,4) = robot.Bodies{i}.CenterOfMass' ;
        Tc(:,:,i) = Tj(:,:,i) * K;

        % orientation of com frame {i} in com frame {i-1}
        if i == 1
           M(:,:,i) = Tc(:,:,i);
        else
           M(:,:,i) = TransInv(Tc(:,:,i-1)) * Tc(:,:,i);
        end

        % screw axes
        w = Tj(1:3,1:3,i) * robot.Bodies{i}.Joint.JointAxis' ; 
        v = cross(-w, Tj(1:3, 4, i));
        S(:,i) = [w;v];
        
        %~~~~~~~~~~~~~~~~~inertial parameters
        switch p.Results.set_type
            case 'point'
                % mass
                mass(i) = p.Results.mass_range(1)*robot.Bodies{i}.Mass;

                % center of mass
                com(:, i) = p.Results.com_range(1)*robot.Bodies{i}.CenterOfMass';

                % inertia
                I{i} = parallel_axis(p.Results.mass_range(1)*robot.Bodies{i}.Inertia, mass(i), com(:,i));  

                % spatial inertias
                G{i} = [   I{i}        zeros(3,3); 
                        zeros(3,3)   mass(i)*eye(3)];
            case 'interval'
                if strcmp(p.Results.add_uncertainty_to, 'all') || (strcmp(p.Results.add_uncertainty_to, 'link') && any(strcmp(p.Results.links_with_uncertainty, robot.Bodies{i}.Name)))
                    % mass
                    lo_mass(i) = p.Results.mass_range(1)*robot.Bodies{i}.Mass;
                    hi_mass(i) = p.Results.mass_range(2)*robot.Bodies{i}.Mass;
                    mass(i) = interval(lo_mass(i), hi_mass(i));

                    % center of mass
                    lo_com(:, i) = p.Results.com_range(1)*robot.Bodies{i}.CenterOfMass';
                    hi_com(:, i) = p.Results.com_range(1)*robot.Bodies{i}.CenterOfMass';

                    % inertia
                    lo_inertia_vec = p.Results.mass_range(1)*robot.Bodies{i}.Inertia;
                    hi_inertia_vec = p.Results.mass_range(2)*robot.Bodies{i}.Inertia;
                    
                    lo_I = parallel_axis(lo_inertia_vec,...
                                         lo_mass(i),... 
                                         lo_com(:, i));
                         
                    hi_I = parallel_axis(hi_inertia_vec,...
                                         hi_mass(i),... 
                                         hi_com(:, i));
                                             
                    % swap high and low values
                    for j = 1:3
                        for k = 1:3
                            a = min(lo_I(j,k), hi_I(j,k));
                            b = max(lo_I(j,k), hi_I(j,k));

                            lo_I(j,k) = a;
                            hi_I(j,k) = b;
                        end
                    end
                    
                    % create interval com
                    a = min(lo_com(:,i), hi_com(:,i));
                    b = max(lo_com(:,i), hi_com(:,i));

                    lo_com(:,i) = a;
                    hi_com(:,i) = b;
                    com(:, i) = interval(lo_com(:, i), hi_com(:, i));

                    % create interval inertia        
                    I{i} = interval(lo_I, hi_I); 

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
            otherwise
                error('Set type not recognized.');
        end
                       
    end
    
    params.set_type = p.Results.set_type;
    params.mass = mass;
    params.com = com;
    params.T0 = T0;
    params.H = H;
    params.P = P;
    params.G = G;
    params.I = I;
    params.M = M;
    params.S = S;
    params.num_bodies = robot.NumBodies;
    params.num_joints = num_joints;
    params.num_q = num_q;
    params.joint_types = joint_types;
    params.q_index = q_index;
%     params.active_joints = active_joints;
%     params.fixed_joints = fixed_joints;
    params.joint_axes = joint_axes;
    params.gravity = [0 0 -9.81];
%     params.use_interval = false;
    params.add_uncertainty_to = p.Results.add_uncertainty_to;
    params.links_with_uncertainty = p.Results.links_with_uncertainty;
    params.mass_range = p.Results.mass_range;
    params.com_range = p.Results.com_range;
    params.zono_order = p.Results.zono_order;
end