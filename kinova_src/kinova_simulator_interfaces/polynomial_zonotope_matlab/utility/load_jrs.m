function [jrs] = load_jrs(q, qd)

    % use initial state to load appropriate jrs
    jrs_path = '/home/ngv-08/code/matlab_code/dynamics-dev/dynamics/rnea/polynomial_zonotope/jrs_saved/'; % where the jrss are stored
    jrs_key_tmp = load([jrs_path, 'c_kvi.mat']);
    jrs_key = jrs_key_tmp.c_kvi;

    time_dim = 6;
    k_dim = 4;
    acc_dim = 4;


    for i = 1:length(qd)
        [~, closest_idx] = min(abs(qd(i) - jrs_key));
        jrs_filename = sprintf('%sJRS_%0.3f.mat', jrs_path, jrs_key(closest_idx));
        jrs_load_tmp = load(jrs_filename);
        jrs_load = jrs_load_tmp.JRS;

        % implementation of lines 4 and 5 of Alg. 2:
        % use A to rotate the cos and sin dimensions of jrs
        % then, slice at correct k^v_i value
        for jrs_idx = 1:length(jrs_load)
            A = [cos(q(i)), -sin(q(i)), 0, 0, 0, 0; sin(q(i)), cos(q(i)), 0, 0, 0, 0; [zeros(4, 2), eye(4)]];
            jrs{jrs_idx, 1}{i, 1} = A*zonotope_slice(jrs_load{jrs_idx}, 5, qd(i));
        end

        %%% HACK?? 
        % alter the second half of the JRSs to account for braking
        % acceleration. assumes ARMTD style traj. parameterization.
        G = jrs{1, 1}{i, 1}.Z(k_dim, 2:end);
        k_idx = find(G~=0);
        if length(k_idx) ~= 1
            error('expected one k-sliceable generator');
        end
        delta_k = G(k_idx);
        c_braking = (0 - qd(i))/0.5; % acceleration to stop from initial speed
        delta_braking = (0 - delta_k)/0.5; % extra acceleration to stop depending on choice of k
        for jrs_idx = 51:length(jrs_load)
            Z = jrs{jrs_idx, 1}{i, 1}.Z;
            Z(acc_dim, 1) = c_braking;
            Z(acc_dim, k_idx + 1) = delta_braking;
            jrs{jrs_idx, 1}{i, 1} = zonotope(Z);
        end
    end
end