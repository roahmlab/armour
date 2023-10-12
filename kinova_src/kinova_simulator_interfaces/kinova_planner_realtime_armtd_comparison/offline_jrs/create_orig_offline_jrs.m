% create_offline_jrs.m
% patrick holmes, 20221025

% purpose: use CORA to precompute JRSs for comparing the new method of
% computing JRSs vs. the old method of using CORA to do the sin/cos stuff

clear; clc;

plot_on = 1;

% define hyperparameters for our particular parameterization:
dim = 6; % 1 (cos), 2 (sin), 3 (q), 4 (dq), 5 (kai), 6 (kvi)
t_plan = 0.5;
t_total = 1;
dt = 0.01;

% generates dynamics parameterized by K
[orig_dyn_zero_to_t_plan, orig_dyn_t_plan_to_t_total] = ...
    generate_orig_parameterized_dynamics(t_plan, t_total);

% described in Sec. 5.5: grid K^v_i space into smaller subintervals,
% compute separate JRSs in each one
n_JRS = 401; % separate initial velocity space (K^v_i) into 401 smaller intervals
c_kvi = linspace(-pi, pi, n_JRS); % centers of initial velocity subintervals
delta_kvi = (c_kvi(2) - c_kvi(1))/2; % subinterval is c_kvi +- delta_kvi
c_kai = 0; % acceleration parameter space (K^a_i) for each JRS centered at 0

% create folder to save precomputed JRSs
if ~exist('../armtd-dev/src/offline_jrs/orig_parameterization', 'dir')
    mkdir('../armtd-dev/src/offline_jrs/orig_parameterization');
end

% save vector of initial velocity subinterval centers
save('../armtd-dev/src/offline_jrs/orig_parameterization/c_kvi.mat', 'c_kvi');

% set options for reachability analysis:
options.timeStep = dt;
options.taylorTerms=5; % number of taylor terms for reachable sets
options.zonotopeOrder= 2; % zonotope order... increase this for more complicated systems.
options.maxError = 1000*ones(dim, 1); % our zonotopes shouldn't be "splitting", so this term doesn't matter for now
options.verbose = false;
% options.uTrans = 0; % we won't be using any inputs, as traj. params specify trajectories
params.U = zonotope([0, 0]);
% options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';
options.alg = 'lin';

for j = 1:n_JRS
    % break JRS computation into two steps...
    tic;
    delta_kai = max(pi/24, abs(c_kvi(j)/3));
    % first, use dyn_zero_to_t_plan dynamics
    params.tStart = 0; % start time
    params.tFinal = t_plan; % end time for these dynamics

    % create initial zonotope for reachability analysis
    % this is described by eqs. (3) and (8)
    x0 = [1; 0; 0; c_kvi(j); c_kai; c_kvi(j)]; % start at q_i = 0, so cos(q_i) = 1 and sin(q_i) = 0
    % use two generators, one for K^a_i and one for K^v_i (eq. 8)
    params.R0 = zonotope([x0, [0; 0; 0; 0; delta_kai; 0], [0; 0; 0; delta_kvi; 0; delta_kvi]]);

    % create system for reachability analysis (1 dummy input)
    % CORA 2018:
    % sys = nonlinearSys(dim, 1, dyn_zero_to_t_plan, options);
    % CORA 2020:
    sys = nonlinearSys(orig_dyn_zero_to_t_plan, dim, 1);


    % compute JRS over t \in [0, t_plan]
    JRS_zero_to_t_plan = reach(sys, params, options);

    % we'll use the last zonotope of JRS_zero_to_t_plan as the initial
    % zonotope for the braking dynamics portion of the JRS. however, we
    % will use the subset of the zonotope corresponding to t = t_plan:
    % (slicing described in Alg. 1)
    % t_plan_slice = zonotope_slice(JRS_zero_to_t_plan{end}{1}, 6, t_plan);
    t_plan_slice = JRS_zero_to_t_plan.timePoint.set{end};
    params.R0 = t_plan_slice;
    params.tStart = t_plan;
    params.tFinal = t_total;

    % create system for reachability analysis (1 dummy input)
    % CORA 2018:
    % sys = nonlinearSys(dim, 1, dyn_t_plan_to_t_total, options);
    % CORA 2020:
    sys = nonlinearSys(orig_dyn_t_plan_to_t_total, dim, 1);

    % compute JRS over t \in [t_plan, t_total]
    JRS_t_plan_to_t_total = reach(sys, params, options);

    % concatenate the full JRS (eq. 9, shown in Fig. 2)
    JRS = [JRS_zero_to_t_plan.timeInterval.set; JRS_t_plan_to_t_total.timeInterval.set];

    % save this JRS
    current_c_kvi = c_kvi(j);
    filename = sprintf('../armtd-dev/src/offline_jrs/orig_parameterization/JRS_%0.3f.mat', current_c_kvi);
    save(filename, 'JRS', 'options', 't_plan', 't_total', 'current_c_kvi');

    % display:
    fprintf('Current initial velocity: %0.3f, ', current_c_kvi);
    toc;

    % plot if you want
    if plot_on
       figure(1); clf; hold on;
       for i = 1:50
           p = plot(JRS{i}, [1, 2], 'b', 'Filled', true);
           p.FaceAlpha = 0.05;
       end
       for i = 51:100
           p = plot(JRS{i}, [1, 2], 'r', 'Filled', true);
           p.FaceAlpha = 0.05;
       end
    end
end
