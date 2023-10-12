function [R, R_t] = get_transforms_from_jrs(jrs, joint_axes)

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
         
        % get zonotope in matrix form
        Z = [jrs{i}.c, jrs{i}.G];
        id_add = jrs{i}.id;
        ngen = size(Z, 2) - 1;
        
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
    end
end
