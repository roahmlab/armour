function [true_params, nominal_params, interval_params, poly_zono_params] = get_all_params(robot_name, use_random)


    % true, nominal, interval, params
    [true_params, nominal_params, interval_params, ~, ~] = get_interval_params(robot_name, use_random);
        
    % intiailize polynomial zonotope params
    poly_zono_params = nominal_params;
    poly_zono_params.use_poly_zono = true;
    
    id = [];
    id_names = [];
    
    for i = 1:poly_zono_params.num_joints  
        % center of mass for each link
        dim = 3;
        inf = interval_params.com.inf(:, i);
        sup = interval_params.com.sup(:, i);
        c = (inf+sup)./2;
        g = zeros(dim, dim);
        for j = 1:dim
           g(j, j) = (sup(j) - inf(j))./2;
        end
        [id, id_names, id_add] = update_id(id, id_names, 3, ['c', num2str(i)]);
        com{i, 1} = polyZonotope_ROAHM(c, g, [], eye(3), id_add);

        % link masses
        % create one generator per interval dimension
        inf = interval_params.mass.inf(i);
        sup = interval_params.mass.sup(i);
        C = (inf + sup)./2 * eye(dim);
        G = [];
        G = cat(3, G, (sup - inf)./2*eye(dim));
        [id, id_names, id_add] = update_id(id, id_names, 1, ['m', num2str(i)]);
        mass{i, 1} = matPolyZonotope_ROAHM(C, G, [], eye(1), id_add);

        % inertia wrt to center of mass frame
        inf = interval_params.I{i}.inf;
        sup = interval_params.I{i}.sup;
        C = (inf + sup)./2;
        G = [];
        for j = 1:dim
            for k = j:dim
                G_tmp = zeros(dim, dim);
                G_tmp(j, k) = (sup(j,k) - inf(j,k))./2;
                if j ~= k
                    G_tmp(k, j) = G_tmp(j, k);
                end
                G = cat(3, G, G_tmp);
            end
        end
        [id, id_names, id_add] = update_id(id, id_names, size(G, 3), ['I', num2str(i)]);
        I{i, 1} = matPolyZonotope_ROAHM(C, G, [], eye(size(G, 3)), id_add);
    end
                  
    poly_zono_params.mass = mass;
    poly_zono_params.com = com;
    poly_zono_params.I = I;
    
    poly_zono_params.id = id;
    poly_zono_params.id_names = id_names;
end
