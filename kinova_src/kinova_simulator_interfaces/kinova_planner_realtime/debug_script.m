% 0. provide proper input data in buffer/armour.in
% 1. run ./compile_debug_script.sh
% 2. run ./test
% comment: you can change the tracking error to 0 for this debug script

close all; clear; clc;

%% initialize robot
robot = importrobot('kinova_without_gripper.urdf');
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', 'all', ...
                           'uncertain_mass_range', [0.97, 1.03]);

A = uarmtd_agent(robot, params,...
                 'animation_set_axes_flag', 0,... 
                 'animation_set_view_flag', 0,...
                 'use_CAD_flag', true,...
                 'add_measurement_noise_', false, ...
                 'measurement_noise_size_', 0);
transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; % matlab doesn't import these from urdf so hard code into class

link_poly_zonotopes = create_pz_bounding_boxes(robot);

%% initialize desired trajectories
% choose random initial conditions and make sure they are aligned with
% the first three rows in buffer/armour.in
q0 = [-1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000 1.0000000000 1.0000000000 1.0000000000]';
qd0 = [1.0000000000 1.0000000000 1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000 -1.0000000000]'; 
qdd0 = [2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000 2.0000000000]';

% choose a random k_range and make sure they are aligned with k_range in
% Parameters.h
k_range = [pi/24, pi/24, pi/24, pi/24, pi/24, pi/24, pi/24]';
% k_range = [pi/48, pi/48, pi/48, pi/48, pi/48, pi/48, pi/60]';

% choose a random point to slice and make sure they are equal to variable
% factors defined in PZ_test.cpp
k = [0.5, 0.6, 0.7, 0.0, -0.5, -0.6, -0.7]';

q1 = q0 + k .* k_range;
qd1 = zeros(7,1);
qdd1 = zeros(7,1);

beta = match_deg5_bernstein_coefficients({q0, qd0, qdd0, q1, qd1, qdd1}, 2);

tspan = linspace(0, 1, 128 + 1);

%% read CUDA output
link_reachset_center = readmatrix('buffer/armour_joint_position_center.out', 'FileType', 'text');
link_reachset_generators = readmatrix('buffer/armour_joint_position_radius.out', 'FileType', 'text');

torque_reachset_center = readmatrix('buffer/armour_constraints.out', 'FileType', 'text');
torque_reachset_radius = readmatrix('buffer/armour_control_input_radius.out', 'FileType', 'text');

%% verification
% figure; view(3); axis equal; hold on; axis on;
% 
% % for tid = 1:128
% tid = 128;
% 
% % choose a random time inside this time interval
% t_lb = tspan(tid);
% t_ub = tspan(tid + 1);
% t = (t_ub - t_lb) * rand + t_lb;
% 
% q = get_desired_traj(beta, t);
% 
% % plot robot
% A.plot_at_time(q);
% 
% % plot link reachsets
% for j = 1:7
%     c = link_reachset_center((tid-1)*7+j, :)';
%     g = link_reachset_generators( ((tid-1)*7+j-1)*3+1 : ((tid-1)*7+j)*3, :);
%     Z = zonotope(c, g);
%     Z_v = vertices(Z)';
%     trisurf(convhulln(Z_v),Z_v(:,1),Z_v(:,2),Z_v(:,3),'FaceColor',[0,0,1],'FaceAlpha',0.1,'EdgeColor',[0,0,1],'EdgeAlpha',0.3);
% end

% end

% [q, qd, qdd] = get_desired_traj(beta, linspace(0,2,128), 2);
% 
% for i = 1:7
%     subplot(3,3,i);
%     plot(linspace(0,1,128), qd(i,:), 'b');
%     hold on;
%     plot(linspace(0,1,128), torque_reachset_center(:,i), 'r');
% end
% sgtitle('desired velocity plot');

figure; hold on;

us = zeros(7,128);
ts = zeros(1,128);
for tid = 1:128
    % choose a random time inside this time interval
    t_lb = tspan(tid);
    t_ub = tspan(tid + 1);
    ts(tid) = (t_ub - t_lb) * rand + t_lb;

    [q, qd, qdd] = get_desired_traj(beta, ts(tid));

    us(:,tid) = rnea(q, qd, qd, qdd, true, params.nominal) + transmision_inertia' .* qdd;
end

u_lb = torque_reachset_center - torque_reachset_radius;
u_ub = torque_reachset_center + torque_reachset_radius;

% there is a better way to do this
for i = 1:7
    subplot(3,3,i);
    hold on;
    plot(ts, us(i,:), 'r');
    plot(ts, u_lb(:,i), 'b');
    plot(ts, u_ub(:,i), 'b');
    title(['link ', num2str(i)]);
    xlabel('time (sec)');
    ylabel('torque (N*m)');
end
sgtitle('sliced torque reachable set');


%% helper functions
function [q, qd, qdd] = get_desired_traj(beta, t, T)
    if nargin < 3
        T = 1;
    end

    [B, dB, ddB] = Bezier_kernel_deg5(t / T);
    
    q = zeros(7,length(t));
    qd = zeros(7,length(t));
    qdd = zeros(7,length(t));
    for j = 1:6
        q = q + beta{j} .* B(:,j)';
        qd = qd + beta{j} .* dB(:,j)';
        qdd = qdd + beta{j} .* ddB(:,j)';
    end

    qd = qd / T;
    qdd = qdd / T / T;
end
