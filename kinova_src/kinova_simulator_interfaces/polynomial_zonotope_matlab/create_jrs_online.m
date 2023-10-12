function [Q_des, Qd_des, Qdd_des, Q, Qd, Qd_a, Qdd_a, R_des, R_t_des, R, R_t, jrs_info] = create_jrs_online(q, dq, ddq, joint_axes, taylor_degree, traj_type, add_ultimate_bound, LLC_info)
% patrick 20220330
% implementing jrs computation that avoids going through CORA"s reachability toolbox
% and relies only on PZ arithmetic.
% can probably be greatly sped up by a lot (all?) computations offline.
% implement for 2 different trajectory parameterizations:
% original ARMTD: q(t) = q_0 + \dot{q}_0*t + 1/2*k*t^2
% Bernstein: q(t) = \sum_i=0^3 \beta_i B_i(t)

% patrick 20221027
% moved file from dynamics_dev to armtd_dev
% continue development on armtd_dev!

%% setup
if ~exist('add_ultimate_bound', 'var')
    add_ultimate_bound = true;
end

make_gens_independent = true; % treat time and error as independent generators
% (we probably don't want to do this until after cos/sin has been taken to get the 
% JRSs as tight as possible.)
if ~make_gens_independent
    disp('Treating time and tracking error as dependent generators. Make sure to remove them before creating constraints.');
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

if ~exist('taylor_degree', 'var')
	taylor_degree = 1; % default uses degree 1 taylor expansion
end

if ~exist('traj_type', 'var')
   traj_type = 'orig'; % use original traj parameterization by default
end

% if ~exist('bernstein_waypoint', 'var')
	bernstein_center = zeros(size(q));
% else
% 	bernstein_center = bernstein_waypoint - q;
% end

if ~exist('LLC_info', 'var')
    % hard coded ultimate bound if not passed in
    ultimate_bound = 0.0191;
    k_r = 10; % assuming K_r = k_r*eye(n_q)
else
	if add_ultimate_bound
	    ultimate_bound = LLC_info.ultimate_bound;
	    k_r = LLC_info.Kr;
    else
        k_r = 0;
	end
end

bernstein_final_range = pi/36*ones(n_q, 1);
% bernstein_final_range = [pi/24; pi/72; pi/24; pi/72; pi/72; pi/72; pi/72];

t_f = 1;
t_p = 0.5;
dt = 0.01;
n_t = t_f/dt;
n_t_p = t_p/t_f*n_t; % last time step of "planning" phase before braking phase

% adjust range that k \in [-1, 1] corresponds to...
n_k = n_q;
c_k_orig = zeros(n_k, 1);
g_k_orig = min(max(pi/24, abs(dq/3)), pi/3); % added critical abs() here!!
% g_k = ones(n_q, 1);

%% initialize
id = [];
id_add = [];
id_names = [];

T = cell(n_t, 1);
K = cell(n_t, 1);
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

%% get PZ versions of k, time, and error:
% get traj param PZs
switch traj_type
    case 'orig'
        for j = 1:n_q
            [id, id_names, id_add] = update_id(id, id_names, 1, ['k', num2str(j)]);
            K{j} = polyZonotope_ROAHM(0, 1, [], 1, id_add);
        end
    case 'bernstein'
        for j = 1:n_q
            [id, id_names, id_add] = update_id(id, id_names, 1, ['k', num2str(j)]);
            % K{j} = polyZonotope_ROAHM(q(j), bernstein_final_range, [], 1, id_add);
            K{j} = polyZonotope_ROAHM(0, 1, [], 1, id_add);
        end
end

% get time subinterval PZs
[id, id_names, id_add] = update_id(id, id_names, n_t, ['t1']);
for i = 1:n_t
	T{i} = polyZonotope_ROAHM(dt*(i-1) + dt/2, dt/2, [], 1, id_add(i, :));
end

% get ultimate bound PZs
if add_ultimate_bound
    for j = 1:n_q
        [id, id_names, id_add] = update_id(id, id_names, 1, ['e1']);
        E_p{j} = polyZonotope_ROAHM(0, ultimate_bound/k_r, [], 1, id_add); % position
        [id, id_names, id_add] = update_id(id, id_names, 1, ['e2']);
        E_v{j} = polyZonotope_ROAHM(0, 2*ultimate_bound, [], 1, id_add); % velocity
    end
else
    for j = 1:n_q
        E_p{j} = 0;
        E_v{j} = 0;
    end
end

%% multiply to get JRS
% some preliminaries for traj computation:
switch traj_type
case 'orig'
	for j = 1:n_q
		Qd_plan{j} = dq(j) + (c_k_orig(j) + g_k_orig(j).*K{j}).*t_p;
		Qdd_brake{j} = (-1).*Qd_plan{j}.*(1/(t_f - t_p));
		Q_plan{j} = q(j) + dq(j).*t_p + 0.5.*(c_k_orig(j) + g_k_orig(j).*K{j}).*t_p^2;
	end
case 'bernstein'
    for j = 1:n_q
       q1{j, 1} = q(j) + bernstein_center(j) + bernstein_final_range(j).*K{j}; % final position is initial position +- k \in [-1, 1]
       dq1 = 0;
       ddq1 = 0;
       beta{j} = match_deg5_bernstein_coefficients({q(j); dq(j); ddq(j); q1{j}; dq1; ddq1});
       alpha{j} = bernstein_to_poly(beta{j}, 5);
    end
end

% main loop:
for i = 1:n_t
	for j = 1:n_q
		switch traj_type
		case 'orig'
			if i <= n_t_p
				Q_des{i}{j, 1} = q(j) + dq(j).*T{i} + 0.5.*(c_k_orig(j) + g_k_orig(j).*K{j}).*T{i}.*T{i};
				Qd_des{i}{j, 1} = dq(j) + (c_k_orig(j) + g_k_orig(j).*K{j}).*T{i};
				Qdd_des{i}{j, 1} = (c_k_orig(j) + g_k_orig(j).*K{j});
			else
				Q_des{i}{j, 1} = Q_plan{j} + Qd_plan{j}.*(T{i} - t_p) + 0.5.*Qdd_brake{j}.*(T{i} - t_p).*(T{i} - t_p);
				Qd_des{i}{j, 1} = Qd_plan{j} + Qdd_brake{j}.*(T{i} - t_p);
				Qdd_des{i}{j, 1} = Qdd_brake{j};
			end
		case 'bernstein'
% 			error('not implemented yet');
            Q_des{i}{j, 1} = 0;
            Qd_des{i}{j, 1} = 0;
            Qdd_des{i}{j, 1} = 0;
            for k = 0:5
                Q_des{i}{j, 1} = Q_des{i}{j, 1} + alpha{j}{k+1}.*T{i}.^k;
                if k > 0
                    Qd_des{i}{j, 1} = Qd_des{i}{j, 1} + k*alpha{j}{k+1}.*T{i}.^(k-1);
                end
                if k > 1
                    Qdd_des{i}{j, 1} = Qdd_des{i}{j, 1} + (k)*(k-1)*alpha{j}{k+1}.*T{i}.^(k-2);
                end
            end
		end 

		% add tracking error
		Q{i}{j, 1} = Q_des{i}{j, 1} + E_p{j};
		Qd{i}{j, 1} = Qd_des{i}{j, 1} + E_v{j};
		Qd_a{i}{j, 1} = Qd_des{i}{j, 1} + k_r.*E_p{j};
		Qdd_a{i}{j, 1} = Qdd_des{i}{j, 1} + k_r.*E_v{j};

		% get rotation matrices:
		[R_des{i}{j, 1}, R_t_des{i}{j, 1}] = get_pz_rotations_from_q(Q_des{i}{j, 1}, joint_axes(:, j), taylor_degree);
		[R{i}{j, 1}, R_t{i}{j, 1}] = get_pz_rotations_from_q(Q{i}{j, 1}, joint_axes(:, j), taylor_degree);
	end
end

%% throw away time/error generators and compress independent generators to just one (if 1D):
if make_gens_independent
	for i = 1:n_t
		for j = 1:n_q
			k_id = id(j);
			Q_des{i}{j, 1} = remove_dependence_and_compress(Q_des{i}{j, 1}, k_id);
			Qd_des{i}{j, 1} = remove_dependence_and_compress(Qd_des{i}{j, 1}, k_id);
			Qdd_des{i}{j, 1} = remove_dependence_and_compress(Qdd_des{i}{j, 1}, k_id);
			Q{i}{j, 1} = remove_dependence_and_compress(Q{i}{j, 1}, k_id);
			Qd{i}{j, 1} = remove_dependence_and_compress(Qd{i}{j, 1}, k_id);
			Qd_a{i}{j, 1} = remove_dependence_and_compress(Qd_a{i}{j, 1}, k_id);
			Qdd_a{i}{j, 1} = remove_dependence_and_compress(Qdd_a{i}{j, 1}, k_id);
			R_des{i}{j, 1} = remove_dependence_and_compress_mat(R_des{i}{j, 1}, k_id);
			R_t_des{i}{j, 1} = remove_dependence_and_compress_mat(R_t_des{i}{j, 1}, k_id);
			R{i}{j, 1} = remove_dependence_and_compress_mat(R{i}{j, 1}, k_id);
			R_t{i}{j, 1} = remove_dependence_and_compress_mat(R_t{i}{j, 1}, k_id);

		end
	end

	% remove all ids except for k
	id(n_q+1:end, :) = [];
	id_names(n_q+1:end, :) = [];
end

% save info
jrs_info.id = id;
jrs_info.id_names = id_names;
jrs_info.k_id = (1:n_q)';
jrs_info.n_t = n_t;
jrs_info.n_q = n_q;
jrs_info.n_k = n_k;
jrs_info.c_k_orig = c_k_orig;
jrs_info.g_k_orig = g_k_orig;
jrs_info.c_k_bernstein = bernstein_center;
jrs_info.g_k_bernstein = bernstein_final_range;

end

function B = remove_dependence_and_compress(A, k_id)
	k_id_idx = (A.id == k_id);
	k_slc_idx = (A.expMat(k_id_idx, :) ~= 0 & all(A.expMat(~k_id_idx, :) == 0, 1)); % should only be one!
	if length(find(k_slc_idx)) > 1
		error('There should only be one fully-k-sliceable generator');
	end
	B = polyZonotope_ROAHM(A.c, A.G(k_slc_idx), sum(abs(A.G(~k_slc_idx))) + sum(abs(A.Grest)), A.expMat(k_id_idx, k_slc_idx), k_id);
	% B = polyZonotope_ROAHM(A.c, A.G(k_slc_idx), [A.G(~k_slc_idx), A.Grest], A.expMat(k_id_idx, k_slc_idx), k_id);
end

function B = remove_dependence_and_compress_mat(A, k_id)
	k_id_idx = (A.id == k_id);
	k_slc_idx = (A.expMat(k_id_idx, :) ~= 0 & all(A.expMat(~k_id_idx, :) == 0, 1)); % should only be one!
% 	if length(find(k_slc_idx)) > 1
% 		error('There should only be one fully-k-sliceable generator');
% 	end
	B = matPolyZonotope_ROAHM(A.C, A.G(:, :, k_slc_idx), cat(3, A.G(:, :, ~k_slc_idx), A.Grest), A.expMat(k_id_idx, k_slc_idx), k_id);
end

