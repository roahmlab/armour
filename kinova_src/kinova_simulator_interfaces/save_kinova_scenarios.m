% visualize and save .csv files for hard kinova scenarios
% you can adjust these scenarios in get_kinova_scenario_info.m

clear; clc;

save_on = true; % if you just want to visualize, set to false.
world_save_dir = './simulator_files/testing/saved_worlds/kinova_hard_scenarios';
if save_on && ~exist(world_save_dir, 'dir')
    mkdir(world_save_dir);
end 

include_base_obstacle = false;
dimension = 3;
verbosity = 10;

% set up a dummy kinova agent for plotting
agent_urdf = 'gen3.urdf';
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
model = create_model_from_urdf(agent_urdf);
model = rmfield(model, 'damping');
model = rmfield(model, 'friction');
params = load_robot_params(robot);
A = uarmtd_agent(robot, model, params);
I = A.get_agent_info ;

% loop through scenarios, plotting and saving
for scenario = 1:7
    clc;
    fprintf("THIS IS SCENARIO %d\n\n", scenario);

    [obstacles, start, goal_radius, goal_type, goal] = get_kinova_scenario_info(scenario);

    W = kinova_world_static('create_random_obstacles_flag', false, 'include_base_obstacle', include_base_obstacle, 'goal_radius', goal_radius, 'N_obstacles',length(obstacles),'dimension',dimension,'workspace_goal_check', 0,...
        'verbose',verbosity, 'start', start, 'goal', goal, 'obstacles', obstacles, 'goal_type', goal_type) ;
    W.setup(I) ;
    W.bounds = [-1 1 -1 1 0 2];

    A.state(A.joint_state_indices) = W.start ;

    % plot and pause;
    figure(1); clf; axis equal ; xlim([-1 1]); ylim([-1 1]); zlim([0 2]); grid on; hold on ; view(3);
    plot(A);
    plot(W);
%     pause;

    % save
    if save_on
        filename = sprintf('%s/scene_%03d.csv', world_save_dir, scenario);
        % create .csv file
        write_fetch_scene_to_csv(W, filename);
    end
end