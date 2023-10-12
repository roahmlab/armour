close all; clear; clc;

%% step 0. generate symbolic version of desired trajectory
% time variable
syms t 'real'

% duration variable
syms T 'real'

% desired initial conditions
syms q0 qd0 qdd0 'real'

% Bezier curve degree
degree = 5;

% Bezier coefficients
b = sym('b', [degree + 1, 1], 'real');

% Bernstein polynomials and their time-derivatives
[B, dB, ddB] = Bezier_kernel_deg5(t / T);

% get trajectories
q = B * b;
qd = dB * b / T;
qdd = ddB * b / T / T;

% compute the actual initial condition
q0_actual = simplify(subs(q, t, 0));
qd0_actual = simplify(subs(qd, t, 0));
qdd0_actual = simplify(subs(qdd, t, 0));

% align the initial conditions
% trying to completely replace b(1) and b(2) but failed
b1_res = solve(q0_actual == q0, b(1));

qd0_actual = subs(qd0_actual, b(1), b1_res);
b2_res = solve(qd0_actual == qd0, b(2));

qdd0_actual = subs(qdd0_actual, [b(1); b(2)], [b1_res; b2_res]);
b3_res = solve(qdd0_actual == qdd0, b(3));

b(1) = b1_res;
b(2) = b2_res;
b(3) = b3_res;

% parameter
syms k 'real'

% parameterize the end position
b(4) = q0 + k;
b(5) = q0 + k;
b(6) = q0 + k;

disp(b);

q = B * b;
qd = dB * b / T;
qdd = ddB * b / T / T;

% scaled variable, guaranteed to range in [0,1]
syms s 'real'
assume(0 <= s); assumeAlso(s <= 1);

% simplify the experssion since qd0 is always with T, qdd0 is always with T^2
syms Tqd0 'real'
syms TTqdd0 'real'

q = subs(q, [t / T; qd0*T; qdd0*T^2], [s; Tqd0; TTqdd0]);
qd = subs(qd, [t / T; qd0*T; qdd0*T^2], [s; Tqd0; TTqdd0]);
qdd = subs(qdd, [t / T; qd0*T; qdd0*T^2], [s; Tqd0; TTqdd0]);

% hence, the following discussion is based on s, but applies to any T > 0

%% step 1. remove variable t in our expressions by replacing it with a fixed (scaled) time inteval [s_lb, s_ub]
syms s_lb s_ub 'real'
assume(s_lb < s_ub);

% [s_lb, s_ub] corresponds to a (scaled time interval) 
% if we divide [0,1] in NUM_TIME_STEPS intervals

% get the coefficient with respect to k
[co, va] = coeffs(q, k);
q_k_dep = co(1);
q_k_indep = co(2);
% so that q = co(1) * k + co(2)

% we first deal with q_k_dep
disp(solve(diff(q_k_dep, s) == 0, s)');
figure; ss = 0:0.01:1;
plot(ss, double(subs(q_k_dep, s, ss)));
xlabel('s'); ylabel('q k dep');
% it can be easily verified that q_k_dep is actually 
% monotonically increasing with t, as a result,
q_k_dep_lb = subs(q_k_dep, s, s_lb);
q_k_dep_ub = subs(q_k_dep, s, s_ub);
q_k_dep_center = 0.5 * (q_k_dep_ub + q_k_dep_lb);
q_k_dep_radius = 0.5 * (q_k_dep_ub - q_k_dep_lb);
% so that q_k_dep is bounded by q_k_dep_center + [-q_k_dep_radius, q_k_dep_radius]

% now let's deal with q_k_indep
q_k_indep_extremas = solve(diff(q_k_indep, s) == 0, s)';
disp(q_k_indep_extremas);
% we can see that q_k_indep has 3 extremas with respect to s
% one of the extrema is just 1 so we can omit that
q_k_indep_extremas = q_k_indep_extremas(3:4);
q_k_indep_extremums = subs(q_k_indep, s, q_k_indep_extremas);

% so the lower bound or the upper bound will only have the following cases:
% (1) q_k_indep(s_lb)
% (2) q_k_indep(s_ub)
% (3) q_k_indep_extremums(1) if s_lb <= q_k_indep_extremas(1) <= s_ub
% (4) q_k_indep_extremums(2) if s_lb <= q_k_indep_extremas(2) <= s_ub
% so all we need to do find out the minimum and maximum of these 4 values
% and treat them as the lower bound and the upper bound, respectively
% we would have to code this logic on our own, let us continue with
% symbolic lower bounds and upper bounds
syms q_k_indep_lb q_k_indep_ub 'real'

q_k_indep_center = 0.5 * (q_k_indep_ub + q_k_indep_lb);
q_k_indep_radius = 0.5 * (q_k_indep_ub - q_k_indep_lb);
% so that q_k_indep is bounded by q_k_indep_center + [-q_k_indep_radius, q_k_indep_radius]

% now we have
% q = q_k_dep * k + q_k_indep
%   = (q_k_dep_center + [-q_k_dep_radius, q_k_dep_radius]) * k + q_k_indep_center + [-q_k_indep_radius, q_k_indep_radius]
%   = q_k_indep_center + q_k_dep_center * k +
%     + ([-q_k_indep_radius, q_k_indep_radius] + [-q_k_dep_radius, q_k_dep_radius] * k)

% For [-q_k_dep_radius, q_k_dep_radius] * k, we have to over-approximate k with k_range
syms k_range 'real'
q_c = q_k_indep_center;
q_G = q_k_dep_center;
q_Grest = q_k_indep_radius + q_k_dep_radius * k_range;
% so that q is bounded by q_c + q_G * k + [-q_Grest, q_Grest]

%% step 2. continue to get cos(q) and sin(q)
% consider the following (first order) Taylor expansion
% cos(q) = cos(q_c + q_G * k + [-q_Grest, q_Grest])
%        = cos(q_c) - sin(c) * (q_G * k + [-q_Grest, q_Grest]) - 
%          - 0.5 * cos(q_c + q_G * k + [-q_Grest, q_Grest]) * (q_G * k + [-q_Grest, q_Grest])^2
%        = cos(q_c) - (sin(c) * q_G) * k + 
%          + (sin(c) * q_G) * [-q_Grest, q_Grest] + 
%          + 0.5 * cos(q_c + q_G * [-k_range, k_range] + [-q_Grest, q_Grest]) * (q_G * [-k_range, k_range] + [-q_Grest, q_Grest])^2

% as a result
cos_q_c = cos(q_c);
cos_q_G = -sin(q_c) * q_G;
syms cos_q_Grest 'real' % cos_q_Grest would have to be computed through interval arithmetic
% so that cos(q) is bounded by cos_q_c + cos_q_G * k + [-cos_q_Grest, cos_q_Grest]

% sin(q) would be similar

%% step 3. repeat the same procedure for qd
% get the coefficient with respect to k
[co, va] = coeffs(qd, k);
qd_k_dep = co(1);
qd_k_indep = co(2);
% so that qd = co(1) * k + co(2)

% we first deal with qd_k_dep
disp(solve(diff(qd_k_dep, s) == 0, s)');
figure; ss = 0:0.01:1;
plot(ss, double(subs(subs(qd_k_dep, T, 2), s, ss)));
xlabel('s'); ylabel('qd k dep');
% it can be easily verified that qd_k_dep is  
% monotonically increasing on [0,  0.5], and
% monotonically decreasing on [0.5,1  ]

% now let's deal with qd_k_indep
qd_k_indep_extremas = solve(diff(qd_k_indep, s) == 0, s)';
disp(qd_k_indep_extremas);
% we can see that q_k_indep has 3 extremas with respect to s
% one of the extrema is just 1 so we can omit that
% the discussion would be similar

%% step 4. repeat the same procedure for qdd
% get the coefficient with respect to k
[co, va] = coeffs(qdd, k);
qdd_k_dep = co(1);
qdd_k_indep = co(2);
% so that qd = co(1) * k + co(2)

% we first deal with qd_k_dep
disp(solve(diff(qdd_k_dep, s) == 0, s)');
figure; ss = 0:0.01:1;
plot(ss, double(subs(subs(qdd_k_dep, T, 2), s, ss)));
xlabel('s'); ylabel('qdd k dep');
% it can be easily verified that qdd_k_dep is  
% monotonically increasing on [0,                0.5 - sqrt(3) / 6], and
% monotonically decreasing on [0.5 - sqrt(3) / 6,0.5 + sqrt(3) / 6], and
% monotonically increasing on [0.5 + sqrt(3) / 6,1                ]

% now let's deal with qdd_k_indep
qdd_k_indep_extremas = solve(diff(qdd_k_indep, s) == 0, s)';
disp(qdd_k_indep_extremas);
% we can see that qdd_k_indep has 2 extremas with respect to s
% the discussion would be similar

