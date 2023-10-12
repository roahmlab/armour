function [q_des, qd_des, qdd_des, q, qd, qd_a, qdd_a, r, c_k, delta_k, id, id_names] = prep_traj(q_0, qd_0, uniform_bound, Kr, id, id_names)

    % use initial state to load appropriate jrs
    jrs_path = '/home/ngv-08/code/matlab_code/dynamics-dev/dynamics/rnea/polynomial_zonotope/jrs_saved/'; % where the jrss are stored
    jrs_key_tmp = load([jrs_path, 'c_kvi.mat']);
    jrs_key = jrs_key_tmp.c_kvi;

    cos_dim = 1;
	sin_dim = 2;
    vel_dim = 3;
    k_dim = 4;
    acc_dim = 4;
    time_dim = 6;

    % number of joints
    num_joints = length(q_0);
    
    % placeholders for acceleration center/generator
    c_k = zeros(num_joints, 1);
    delta_k = zeros(num_joints, 1);
        
    %% load joint reachable sets
    for i = 1:length(qd_0)
        % load closest jrs
        [~, closest_idx] = min(abs(qd_0(i) - jrs_key));
        jrs_filename = sprintf('%sJRS_%0.3f.mat', jrs_path, jrs_key(closest_idx));
        jrs_load_tmp = load(jrs_filename);
        jrs_load = jrs_load_tmp.JRS;

        % Alg. 2: lines 4-5 (desired trajectory)
        for jrs_idx = 1:length(jrs_load)
            % create rotation matrix from cos and sin dimensions of jrs
            A = [cos(q_0(i)), -sin(q_0(i)), 0, 0, 0, 0; sin(q_0(i)), cos(q_0(i)), 0, 0, 0, 0; [zeros(4, 2), eye(4)]];
            
            % slice at correct k^v_i value, then rotate using A
            jrs{jrs_idx, 1}{i, 1} = A*zonotope_slice(jrs_load{jrs_idx}, 5, qd_0(i));
        end

        % save interval information
        G = jrs{1, 1}{i, 1}.Z(k_dim, 2:end);
        k_idx = find(G~=0);
                     
        if length(k_idx) ~= 1
            error('expected one k-sliceable generator');
        end
        
        c_k(i,1) = jrs{1,1}{i,1}.Z(k_dim, 1);
        delta_k(i,1) = G(k_idx);
        
        %%% HACK: alter the second half of the JRSs to account for braking
        % acceleration. assumes ARMTD style traj. parameterization.
        
        % acceleration to stop from initial speed
        c_braking = (0 - qd_0(i))/0.5; 
        
        % extra acceleration to stop depending on choice of k
        delta_braking = (0 - delta_k(i,1))/0.5; 
   
        for jrs_idx = 51:length(jrs_load)
            Z = jrs{jrs_idx, 1}{i, 1}.Z;
            Z(acc_dim, 1) = c_braking;
            Z(acc_dim, k_idx + 1) = delta_braking;
            jrs{jrs_idx, 1}{i, 1} = zonotope(Z);    
        end     
    end
        
    %% create error zonotopes
    % controller gainz, brah
    gain = Kr(1,1);
    
    % max position error
    e = uniform_bound / gain;
    
    % max velocity erroy
    d = 2*uniform_bound;
    
    % position error generators
    G_pos_err = zeros(num_joints, 2);
    G_pos_err(1,1) = (cos(0) - cos(e))/2;
    G_pos_err(2,2) = (sin(e) - sin(-e))/2;

    % velocity error generator
    G_vel_err = d;    
    
    % label position error coefficients as 'e'
    [~, pos_err_ngens] = size(G_pos_err);
    [id, id_names, pos_err_id_add] = update_id(id, id_names, pos_err_ngens, ['e', num2str(1)]);
    
    % label velocity error coefficients as 'd'
    [~, vel_err_ngens] = size(G_vel_err);
    [id, id_names, vel_err_id_add] = update_id(id, id_names, vel_err_ngens, ['d', num2str(1)]);    
 
    % label gain*(position error) coefficients as 'p'
    [id, id_names, gain_pos_err_id_add] = update_id(id, id_names, 1, ['p', num2str(1)]);
 
    % label gain*(velocity error) coefficients as 'v'
    [id, id_names, gain_vel_err_id_add] = update_id(id, id_names, 1, ['v', num2str(1)]);
    
    % modified error labels
    r_id_add = [vel_err_id_add; gain_pos_err_id_add];    
    r_ngens = length(r_id_add);
    
    %% create trajectories
    for i = 1:length(jrs) % for each time step
        for j = 1:num_joints
            % get zonotope as matrix
            Z = jrs{i,1}{j}.Z;
            
            % remove zero generators
            C = Z(:,1);
            G = Z(:, 2:end);
            G(:, ~any(G)) = [];
            
            Z = [C, G];
            
            % number of generators
            ngen = size(Z, 2) - 1;
                        
            % find 'k' dependent generator
            k_idx = find(Z(k_dim, 2:end) ~= 0);
            if length(k_idx) == 0
                warning('no k-dependence detected');
                [id, id_names, id_add] = update_id(id, id_names, ngen, ['j', num2str(j)]);
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
                [id, id_names, id_add1] = update_id(id, id_names, 1, ['k', num2str(j)]);
                
                % label the rest of the coefficients as 'j' for JRS
                [id, id_names, id_add2] = update_id(id, id_names, ngen-1, ['j', num2str(j)]);
                id_add = [id_add1; id_add2];                                
            end

            %% desired trajectories
            % save desired positions
            q_des{i, 1}{j,1} = polyZonotope_ROAHM(Z(:, 1), Z(:, 2:end), [], eye(ngen), id_add); % save whole trajectory for slicing later
            
            %%% doing it this way might be unneccessary, but simplifies the
            %%% number of ids for the early zonotopes that we encounter
            
            % save desired velocities
            vel_C = Z(vel_dim, 1);
            vel_G = Z(vel_dim, 2:end);
            
            vel_idx = find(vel_G);
            vel_G = vel_G(:, vel_idx);
            [~, vel_ngen] = size(vel_G);
            qd_des{i, 1}{j,1} = polyZonotope_ROAHM(vel_C, vel_G, [], eye(vel_ngen), id_add(vel_idx'));
            
            % save desired accelerations
            acc_C = Z(acc_dim, 1);
            acc_G = Z(acc_dim, 2:end);
            
            acc_idx = find(acc_G);
            acc_G = acc_G(:, acc_idx);
            [~, acc_ngen] = size(acc_G);
            qdd_des{i, 1}{j,1} = polyZonotope_ROAHM(acc_C, acc_G, [], eye(acc_ngen), id_add(acc_idx')); 
            
            %% actual trajectories
            q{i, 1}{j,1} = polyZonotope_ROAHM(Z(:, 1), [Z(:, 2:end), G_pos_err], [], eye(ngen + pos_err_ngens), [id_add; pos_err_id_add]); % save whole trajectory for slicing later
            qd{i, 1}{j,1} = polyZonotope_ROAHM(vel_C, [vel_G, G_vel_err], [], eye(vel_ngen + vel_err_ngens), [id_add(vel_idx'); vel_err_id_add]);
            
            %% modified trajectories
            qd_a{i, 1}{j,1} = polyZonotope_ROAHM(vel_C, [vel_G, gain*e], [], eye(vel_ngen+1), [id_add(vel_idx'); gain_pos_err_id_add]);
            qdd_a{i, 1}{j,1} = polyZonotope_ROAHM(acc_C, [acc_G, gain*d], [], eye(acc_ngen+1), [id_add(acc_idx'); gain_vel_err_id_add]);
            r{i, 1}{j,1} = polyZonotope_ROAHM(0, [d, gain*e], [], eye(r_ngens), r_id_add);
        end
    end
end

