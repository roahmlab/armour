%% description
% This script generates random worlds to be used for evaluating aRmTD
% against other motion planners
%
% Authors: Patrick Holmes
% Created: 11 November 2019

clear ; clc ; figure(1); clf; view(3); grid on;

%% user parameters
world_save_dir = 'saved_worlds/random';
if ~exist(world_save_dir, 'dir')
    mkdir(world_save_dir);
end

N_obstacle_min = 13 ;
N_obstacle_max = 40 ;
N_obstacle_delta = 3 ;
N_worlds_per_obstacle = 10;
dimension = 3 ;
nLinks = 3 ;
verbosity = 10 ;
allow_replan_errors = true ;
t_plan = 0.5 ;
time_discretization = 0.01 ;
T = 1 ;
use_cuda_flag = false;
agent_move_mode = 'direct' ; % pick 'direct' or 'integrator'

% kinova
robot = importrobot('gen3.urdf');
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
model = create_model_from_urdf('gen3.urdf');
model = rmfield(model, 'transmissionInertia');
model = rmfield(model, 'friction');
model = rmfield(model, 'damping');
params = load_robot_params(robot);
joint_speed_limits = [-1.3963, -1.3963, -1.3963, -1.3963, -1.2218, -1.2218, -1.2218;
                       1.3963,  1.3963,  1.3963,  1.3963,  1.2218,  1.2218,  1.2218]; % matlab doesn't import these from urdf
joint_input_limits = [-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4;
                       56.7,  56.7,  56.7,  56.7,  29.4,  29.4,  29.4]; % matlab doesn't import these from urdf
A = uarmtd_agent(robot, model, params, ...
                     'move_mode', agent_move_mode,...
                     'joint_speed_limits', joint_speed_limits, ...
                     'joint_input_limits', joint_input_limits);
A.LLC = uarmtd_robust_CBF_LLC();

%% automated from here

for i = N_obstacle_min:N_obstacle_delta:N_obstacle_max
    for j = 1:N_worlds_per_obstacle
    
        % use this to start from random start config:
        W = kinova_world_static('include_base_obstacle', 1, 'goal_radius', pi/30, 'N_random_obstacles',i,'dimension',dimension,'workspace_goal_check', 0,...
            'verbose',verbosity, 'creation_buffer', 0.075, 'base_creation_buffer', 0.075) ;

        % set up world using arm
        I = A.get_agent_info ;
        W.setup(I)

        % place arm at starting configuration
        A.state(A.joint_state_indices) = W.start ;
        
        filename = sprintf('%s/scene_%03d_%03d.csv', world_save_dir, i, j);

        % create .csv file
        write_fetch_scene_to_csv(W, filename);

    end
end