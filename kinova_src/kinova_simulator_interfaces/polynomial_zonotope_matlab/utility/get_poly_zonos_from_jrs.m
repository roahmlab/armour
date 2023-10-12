function [R, R_t, joint_vel, joint_acc, id, id_names] = get_poly_zonos_from_jrs(jrs, joint_axes, id, id_names)

	cos_dim = 1;
	sin_dim = 2;
	vel_dim = 3;
	acc_dim = 4;
	k_dim = 4;

    for i = 1:length(jrs) % for each joint
    	rotation_axis = joint_axes(:, i);
    	% use axis-angle formula to get rotation matrix:
    	% https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    	e = rotation_axis./norm(rotation_axis);
    	U = [0 -e(3) e(2);...
    	     e(3) 0 -e(1);...
    	     -e(2) e(1) 0];
    	Z = jrs{i}.Z;
        ngen = size(Z, 2) - 1;
        % find 'k' dependent generator
    	k_idx = find(Z(k_dim, 2:end) ~= 0);
        if length(k_idx) == 0
            warning('no k-dependence detected');
            [id, id_names, id_add] = update_id(id, id_names, ngen, ['j', num2str(i)]);
        elseif length(k_idx) ~= 1
    		error('should only be one non-zero generator in k-dimension');
        else
           % make 'k' dependent generator the first one.
            if k_idx ~= 1
                tmp = Z(:, 2);
                Z(:, 2) = Z(:, k_idx+1);
                Z(:, k_idx+1) = tmp;
            end
            % label the first generator coefficient as 'k'-dependent
            [id, id_names, id_add1] = update_id(id, id_names, 1, ['k', num2str(i)]);
            % label the rest of the coefficients as 'j' for JRS
            [id, id_names, id_add2] = update_id(id, id_names, ngen-1, ['j', num2str(i)]);
            id_add = [id_add1; id_add2];
        end
        G = [];
        G_t = [];
        % create rotation matrices from cos and sin dimensions
    	for j = 1:size(Z, 2)
    		cq = Z(cos_dim, j);
    		sq = Z(sin_dim, j);
    		if j == 1
    			C = eye(3) + sq*U + (1 - cq)*U^2;
    			C_t = C';
            else
                % for generators, create 3x3xn array, where n is number of
                % generators in the jrs
                G_tmp = sq*U + -cq*U^2;
                G = cat(3, G, G_tmp);
                G_t = cat(3, G_t, G_tmp');
			end
        end

    	R{i, 1} = matPolyZonotope_ROAHM(C, G, [], eye(ngen), id_add);
    	R_t{i, 1} = matPolyZonotope_ROAHM(C_t, G_t, [], eye(ngen), id_add);
        joint_vel{i, 1} = polyZonotope_ROAHM(Z(vel_dim, 1), Z(vel_dim, 2:end), [], eye(ngen), id_add);
        joint_acc{i, 1} = polyZonotope_ROAHM(Z(acc_dim, 1), Z(acc_dim, 2:end), [], eye(ngen), id_add);
    end
end