function [orig_dyn_zero_to_t_plan, orig_dyn_t_plan_to_t_total] = generate_orig_parameterized_dynamics(t_plan,t_total)
%GENERATE_ORIG_PARAMETERIZED_DYNAMICS prepares dynamics for CORA 2020 with the
%particular parameterization described in Section 2 of the (original) paper.
%   [dyn_zero_to_t_plan, dyn_t_plan_to_t_total] = generate_parameterized_dynamics(t_plan, t_total)
%   starting from initial velocity k_i^v, we use constant acceleration over
%   t \in [0, t_plan]. 
%   then, we define trajectories with a failsafe
%   (braking) maneuver from the peak speed over t \in [t_plan, t_total].

if ~exist('../armtd-dev/src/offline_jrs/orig_parameterization/jrs_dynamics', 'dir')
   mkdir('../armtd-dev/src/offline_jrs/orig_parameterization/jrs_dynamics')
end

syms cqi sqi qi dqi kai kvi real;
syms udummy real; % CORA will require these arguments, but we won't use them.
x = [cqi; sqi; qi; dqi; kai; kvi];

% these dynamics are written in eqs. (2) and (5) in the paper
ddqi = kai;
dcqi = -sqi*dqi;
dsqi = cqi*dqi;
% dqi = q_i_dot;
dkai = 0;
dkvi = 0;

dx = [dcqi; dsqi; dqi; ddqi; dkai; dkvi];
orig_dyn_zero_to_t_plan = matlabFunction(dx, 'File', '../armtd-dev/src/offline_jrs/orig_parameterization/jrs_dynamics/orig_dyn_zero_to_t_plan', 'vars', {x, udummy});

% now we specify braking dynamics on t \in [t_plan, t_total]
t_to_stop = t_total - t_plan;
q_i_dot_pk = kvi + kai*t_plan;
braking_acceleration = (0 - q_i_dot_pk)/t_to_stop; % brake to 0 velocity from q_i_dot_pk in t_to_stop seconds
ddqi = braking_acceleration;
% q_i_dot = q_i_dot_pk + braking_acceleration*(t - t_plan);
dcqi = -sqi*dqi;
dsqi = cqi*dqi;
% dqi = q_i_dot;

dx = [dcqi; dsqi; dqi; ddqi; dkai; dkvi];
orig_dyn_t_plan_to_t_total = matlabFunction(dx, 'File', '../armtd-dev/src/offline_jrs/orig_parameterization/jrs_dynamics/orig_dyn_t_plan_to_t_total', 'vars', {x, udummy});

end
