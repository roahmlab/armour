%% description
% This script defines a simple world with user defined obstacles and runs
% the ARMOUR planner on them. It then saves information on how well each the
% planner performed in each trial.
%
% Authors: Bohao Zhang (adapted from Patrick Holmes code)
% Created 25 November 2022

initialize_script_path = matlab.desktop.editor.getActiveFilename;
cd(initialize_script_path(1:end-28));

close all; clear; clc;

%% user parameters
goal_type = 'configuration'; % pick 'end_effector_location' or 'configuration'
goal_radius = pi/30;
dimension = 3 ;
verbosity = 10;

%%% for planner
traj_type = 'bernstein'; % pick 'orig' (ARMTD) or 'bernstein' (ARMOUR)
use_cuda_flag = true;

%%% for agent
agent_urdf = 'kinova_with_dumbbell.urdf';

add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.
uncertain_mass_range = [0.97, 1.03];

agent_move_mode = 'direct' ; % pick 'direct' or 'integrator'
use_CAD_flag = false; % plot robot with CAD or bounding boxes

%%% for LLC
use_robust_input = true;
LLC_V_max = 5e-5;

%%% for HLP
if_use_RRT = true;
HLP_grow_tree_mode = 'new' ; % pick 'new' or 'keep'
plot_waypoint_flag = true ;
plot_waypoint_arm_flag  = true ;
lookahead_distance = 0.1 ;

% plotting
plot_while_running = true ;

% simulation
max_sim_time = 172800 ; % 48 hours
max_sim_iter = 600 ;
stop_threshold = 5 ; % number of failed iterations before exiting

%%% for world
start = [0.28159142, 0.6, 3.6350845 , 4.73481646 - 2*pi, 2.55565072, 0.83405794, 2.05711487]'; % start configuration
goal = [4.56258737 - 2*pi, 0.4, 3.39192523, 4.20470761 - 2*pi, 2.86698255, 0.41219441, 6.05075981 - 2*pi]'; % goal configuration
obstacles{1} = box_obstacle_zonotope('center', [0.43014 + 0.05;  0.16598;  0.172],...
                                     'side_lengths', [0.33; 0.33; 0.33]) ;
obstacles{2} = box_obstacle_zonotope('center', [-0.3;  0;  0.5],...
                                     'side_lengths', [0.01; 2; 2]) ;
obstacles{3} = box_obstacle_zonotope('center', [0.8; 0; 1],...
                                     'side_lengths', [2; 2; 0.01]) ;
obstacles{4} = box_obstacle_zonotope('center', [1; -0.6; 0],...
                                     'side_lengths', [2; 0.01; 2]) ;
obstacles{5} = box_obstacle_zonotope('center', [1; 0.81; 0],...
                                     'side_lengths', [2; 0.01; 1]) ;
%% robot params:
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot, ...
                           'add_uncertainty_to', add_uncertainty_to, ...
                           'links_with_uncertainty', links_with_uncertainty,...
                           'uncertain_mass_range', uncertain_mass_range);
joint_speed_limits = [-1.3963, -1.3963, -1.3963, -1.3963, -1.2218, -1.2218, -1.2218;
                       1.3963,  1.3963,  1.3963,  1.3963,  1.2218,  1.2218,  1.2218]; % matlab doesn't import these from urdf so hard code into class
joint_input_limits = [-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4;
                       56.7,  56.7,  56.7,  56.7,  29.4,  29.4,  29.4]; % matlab doesn't import these from urdf so hard code into class
transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; % matlab doesn't import these from urdf so hard code into class
M_min_eigenvalue = 5.095620491878957; % matlab doesn't import these from urdf so hard code into class

%% automated from here
if plot_while_running
    figure(1); clf; view(3); grid on;
end

% run loop
tic;
W = kinova_world_static('create_random_obstacles_flag', false, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
                        'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;

% create arm agent
A = uarmtd_agent(robot, params,...
                 'verbose', verbosity,...
                 'animation_set_axes_flag', 0,... 
                 'animation_set_view_flag', 0,...
                 'move_mode', agent_move_mode,...
                 'use_CAD_flag', use_CAD_flag,...
                 'joint_speed_limits', joint_speed_limits, ...
                 'joint_input_limits', joint_input_limits, ...
                 'add_measurement_noise_', false, ...
                 'measurement_noise_size_', 0,...
                 'M_min_eigenvalue', M_min_eigenvalue, ...
                 'transmision_inertia', transmision_inertia);

% LLC
if use_robust_input
    A.LLC = uarmtd_robust_CBF_LLC('verbose', verbosity, ...
                                  'use_true_params_for_robust', false, ...
                                  'V_max', LLC_V_max, ...
                                  'if_use_mex_controller', true);
else
    A.LLC = uarmtd_nominal_passivity_LLC('verbose', verbosity);
end

A.LLC.setup(A);

P = uarmtd_planner('verbose', verbosity, ...
                   'first_iter_pause_flag', false, ...
                   'use_q_plan_for_cost', true, ...
                   'input_constraints_flag', true, ...
                   'use_robust_input', use_robust_input, ...
                   'traj_type', traj_type, ...
                   'use_cuda', use_cuda_flag, ...
                   'plot_HLP_flag', true) ;

if if_use_RRT
    P.HLP = arm_end_effector_RRT_star_HLP('plot_waypoint_flag',plot_waypoint_flag,...
                                          'plot_waypoint_arm_flag',plot_waypoint_arm_flag,...
                                          'grow_tree_mode',HLP_grow_tree_mode,...
                                          'buffer',0.25, ...
                                          'grow_tree_mode', 'seed') ;
end

% set up world using arm
I = A.get_agent_info ;
W.setup(I) ;
W.bounds = [-1 1 -1 1 0 2];

% place arm at starting configuration
A.state(A.joint_state_indices) = W.start ;

% create simulator
S = simulator_armtd(A,W,P, ...
                    'verbose', verbosity, ...
                    'stop_threshold', stop_threshold, ...
                    'plot_while_running', plot_while_running,...
                    'allow_replan_errors',true,...
                    'max_sim_time',max_sim_time,...
                    'max_sim_iterations',max_sim_iter,...
                    'stop_sim_when_ultimate_bound_exceeded', use_robust_input) ; 

% %% plotting
if plot_while_running
    figure(1) ; clf ; axis equal ; xlim([-1 1]); ylim([-1 1]); zlim([0 2]); grid on; hold on ;

    if dimension == 3
        view(3);
    end
    
    plot(A);
    plot(W);
end

% run simulation
summary = S.run();

