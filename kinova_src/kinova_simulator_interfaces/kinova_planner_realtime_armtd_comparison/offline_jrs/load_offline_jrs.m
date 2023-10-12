function [Q_des, Qd_des, Qdd_des, Q, Qd, Qd_a, Qdd_a, R_des, R_t_des, R, R_t, jrs_info] = load_offline_jrs(q, dq, ddq, joint_axes, add_ultimate_bound)
% patrick 20221025
% function to load offline JRSs as was originally done for ARMTD
% modeled after "create_jrs_online" which replaced the offline jrs loading
% original ARMTD parameterization: q(t) = q_0 + \dot{q}_0*t + 1/2*k*t^2

% use initial state to load appropriate jrs
jrs_path = '../armtd-dev/src/offline_jrs/orig_parameterization/'; % where the jrss are stored
jrs_key_tmp = load([jrs_path, 'c_kvi.mat']);
jrs_key = jrs_key_tmp.c_kvi;

cos_dim = 1;
sin_dim = 2;
pos_dim = 3;
vel_dim = 4;
k_a_dim = 5;
k_v_dim = 6;

%% setup
id = [];
id_names = [];
% hard coded ultimate bound... update to be inputs:
ultimate_bound = 0.0191;
k_r = 10; % assuming K_r = k_r*eye(n_q)
if ~exist('add_ultimate_bound', 'var')
    add_ultimate_bound = false;
end
if add_ultimate_bound
	E_p = polyZonotope_ROAHM(0, [], ultimate_bound/k_r);
	E_v = polyZonotope_ROAHM(0, [], 2*ultimate_bound);
	sin_E_p_g = sin(ultimate_bound/k_r);
	cos_E_p_g = sin_E_p_g;
else
	E_p = polyZonotope_ROAHM(0);
	E_v = polyZonotope_ROAHM(0);
	sin_E_p_g = 0;
	cos_E_p_g = 0;
end

q = q(:);
dq = dq(:);
n_q = length(q);

if ~exist('ddq', 'var')
    ddq = zeros(n_q, 1);
end

if ~exist('joint_axes', 'var')
	joint_axes = [zeros(2, length(q)); ones(1, length(q))]; % default assumes all rotations about z axis
end

traj_type = 'orig'; % use original traj parameterization by default

t_f = 1;
t_p = 0.5;
dt = 0.01;
n_t = t_f/dt;
n_t_p = t_p/t_f*n_t; % last time step of "planning" phase before braking phase

% adjust range that k \in [-1, 1] corresponds to...
n_k = n_q;
% c_k = zeros(n_k, 1);
% g_k = min(max(pi/24, dq/3), pi/3);
% g_k = ones(n_q, 1);

%% initialize
id = [];
id_add = [];
id_names = [];

Q_des = cell(n_t, 1);
Qd_des = cell(n_t, 1);
Qdd_des = cell(n_t, 1);
Q = cell(n_t, 1);
Qd = cell(n_t, 1);
Qd_a = cell(n_t, 1);
Qdd_a = cell(n_t, 1);
R_des = cell(n_t, 1);
R_t_des = cell(n_t, 1);
R = cell(n_t, 1);
R_t = cell(n_t, 1);

%% load joint reachable sets
for j = 1:length(dq)
    % load closest jrs
    [~, closest_idx] = min(abs(dq(j) - jrs_key));
    jrs_filename = sprintf('%sJRS_%0.3f.mat', jrs_path, jrs_key(closest_idx));
    jrs_load_tmp = load(jrs_filename);
    jrs_load = jrs_load_tmp.JRS;

    % Alg. 2: lines 4-5 (desired trajectory)
    for i = 1:length(jrs_load)
        % create rotation matrix from cos and sin dimensions of jrs
        A = [cos(q(j)), -sin(q(j)), 0, 0, 0, 0; sin(q(j)), cos(q(j)), 0, 0, 0, 0; [zeros(4, 2), eye(4)]];
        
        % slice at correct k^v_i value, then rotate using A
        jrs{i, 1} = A*zonotope_slice(jrs_load{i}, 6, dq(j));

        % shift position center
        jrs{i, 1} = jrs{i, 1} + [0; 0; q(j); 0; 0; 0];
    end

    % save k interval information
    G = jrs{1, 1}.Z(k_a_dim, 2:end);
    k_idx = find(G~=0);
                 
    if length(k_idx) ~= 1
        error('expected one k-sliceable generator');
    else
    	% label the first generator coefficient as 'k'-dependent
    	[id, id_names, id_add] = update_id(id, id_names, 1, ['k', num2str(j)]);
    end
    
    c_k(j,1) = jrs{1,1}.Z(k_a_dim, 1);
    g_k(j,1) = G(k_idx);

    % create poly zonos from JRS
    for i = 1:length(jrs_load)
    	% get zonotope as matrix
    	Z = jrs{i,1}.Z;
    	
    	% remove zero generators
    	C = Z(:,1);
    	G = Z(:, 2:end);
    	G(:, ~any(G)) = [];
    	
    	Z = [C, G];
    	
    	% number of generators
    	ngen = size(Z, 2) - 1;
    	            
    	% find 'k' dependent generator
    	k_idx = find(Z(k_a_dim, 2:end) ~= 0);
    	if length(k_idx) == 0
    	    error('no k-dependence detected');
    	elseif length(k_idx) ~= 1
    	    error('should only be one non-zero generator in k-dimension');
    	else
    	   % make 'k' dependent generator the first one.
    	    if k_idx ~= 1
    	        tmp = Z(:, 2);
    	        Z(:, 2) = Z(:, k_idx+1);
    	        Z(:, k_idx+1) = tmp;
    	    end                            
    	end

    	%% desired trajectories
    	% save desired positions
    	Q_des{i, 1}{j, 1} = polyZonotope_ROAHM(Z(pos_dim, 1), Z(pos_dim, 2), Z(pos_dim, 3:end), 1, id_add); % save whole trajectory for slicing later
    	Qd_des{i, 1}{j, 1} = polyZonotope_ROAHM(Z(vel_dim, 1), Z(vel_dim, 2), Z(vel_dim, 3:end), 1, id_add); % save whole trajectory for slicing later
    	% account for braking in JRS trajectories:
    	if i <= 50
	    	Qdd_des{i, 1}{j, 1} = polyZonotope_ROAHM(Z(k_a_dim, 1), Z(k_a_dim, 2), [], 1, id_add);
	    else
		    % acceleration to stop from initial speed
		    c_braking = (0 - dq(j))/0.5; 
    
		    % extra acceleration to stop depending on choice of k
		    delta_braking = (0 - (0.5*g_k(j,1)))/0.5; 
	    	Qdd_des{i, 1}{j, 1} = polyZonotope_ROAHM(c_braking, delta_braking, [], 1, id_add);
    	end

    	%% add tracking error
    	Q{i, 1}{j, 1} = Q_des{i, 1}{j, 1} + E_p;
    	Qd{i, 1}{j, 1} = Qd_des{i}{j, 1} + E_v;
		Qd_a{i, 1}{j, 1} = Qd_des{i}{j, 1} + k_r.*E_p;
		Qdd_a{i, 1}{j, 1} = Qdd_des{i}{j, 1} + k_r.*E_v;
    	
    	%% desired rotations
    	rotation_axis = joint_axes(:, j);
    	% use axis-angle formula to get rotation matrix:
    	% https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    	e = rotation_axis./norm(rotation_axis);
    	U = [0 -e(3) e(2);...
    	     e(3) 0 -e(1);...
    	     -e(2) e(1) 0];
	    R_G = [];
        R_G_t = [];
        % create rotation matrices from cos and sin dimensions
    	for l = 1:size(Z, 2)
    		cq = Z(cos_dim, l);
    		sq = Z(sin_dim, l);
    		if l == 1
    			R_C = eye(3) + sq*U + (1 - cq)*U^2;
    			R_C_t = R_C';
            else
                % for generators, create 3x3xn array, where n is number of
                % generators in the jrs
                G_tmp = sq*U + -cq*U^2;
                R_G = cat(3, R_G, G_tmp);
                R_G_t = cat(3, R_G_t, G_tmp');
			end
        end
        R_des{i, 1}{j, 1} = matPolyZonotope_ROAHM(R_C, R_G(:, :, 1), R_G(:, :, 2:end), 1, id_add);
    	R_t_des{i, 1}{j, 1} = matPolyZonotope_ROAHM(R_C_t, R_G_t(:, :, 1), R_G_t(:, :, 2:end), 1, id_add);

    	% add ultimate bound to rotation matrices
    	R_G = cat(3, R_G, sin_E_p_g*U);
    	R_G = cat(3, R_G, -cos_E_p_g*U^2);
    	R_G_t = cat(3, R_G_t, (sin_E_p_g*U)');
    	R_G_t = cat(3, R_G_t, (-cos_E_p_g*U^2)');
        R{i, 1}{j, 1} = matPolyZonotope_ROAHM(R_C, R_G(:, :, 1), R_G(:, :, 2:end), 1, id_add);
    	R_t{i, 1}{j, 1} = matPolyZonotope_ROAHM(R_C_t, R_G_t(:, :, 1), R_G_t(:, :, 2:end), 1, id_add);
    end  
end

jrs_info.id = id;
jrs_info.id_names = id_names;
jrs_info.k_id = (1:n_q)';
jrs_info.n_t = n_t;
jrs_info.n_q = n_q;
jrs_info.n_k = n_k;
jrs_info.c_k = c_k;
jrs_info.g_k = g_k;


end