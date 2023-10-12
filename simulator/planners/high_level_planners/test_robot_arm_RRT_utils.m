%% Basic RRT functionality tests
% Daphna Raz 2019

%% 2D tests
% dimension implicity set by size of goal
goal = [0 pi/2]';
RRT_object = robot_arm_RRT_HLP_2D('goal', goal);

%% test generate_random_config(HLP, bounds)
bounds = [0 -2*pi;
          pi 2*pi;];

q_rand = RRT_object.generate_random_config(bounds);

for i = 1:RRT_object.dimension
    assert(q_rand(i) > bounds(1,i) & q_rand(i) < bounds(2,i))
end
    

%% test find_nearest_node(HLP, q_rand)
node1 = [0 0]';
node2 = [0 pi/2]';
node3 = [2*pi 2*pi]';

nodes = [node1 node2 node3];

RRT_object.nodes = nodes;
RRT_object.node_parent_idx = [0 1 2];
q_near_idx = RRT_object.find_nearest_node(goal);
assert(q_near_idx == 2);

q_near = RRT_object.nodes(1:RRT_object.dimension, q_near_idx);

%% test generate_new_config(HLP, q_rand, q_near_idx)

q_new = RRT_object.generate_new_config(q_rand,q_near_idx);

% verify that q_rand, q_near, and q_new are collinear
C = [q_rand-q_near q_near-q_new];

assert(rank(C) < 2);

%% Test collision checking 

%% user parameters
% configuration for the 2-link, 2-DOF, 2-D arm
q = rand(2,1) ;

% obstacle
obstacle_size = 0.3 ; % this is a scalar value
obstacle_offset = [0.5 ; 0.5] ;


%% automated from here
% create arm
%
clear ; clc ;

goal = [0 pi/2]';
RRT_object = robot_arm_RRT_HLP_2D('goal', goal);
%% user parameters
% number of obstacles
N_obstacles = 5 ;

% dimension of world
dimension = 2 ; % 2 or 3

% world verbosity
world_verbosity = 10 ;

A = robot_arm_2D_2DOF

%% automated from here

% make world
W = arm_world_static('N_obstacles',N_obstacles,'dimension',dimension,...
    'verbose',world_verbosity) ;

% set up world using arm
I = A.get_agent_info ;
W.setup(I)

% place arm at starting configuration
A.state(A.joint_state_indices) = W.start ;

% TODO get current state 
q_check = W.start

O = W.obstacles 

out = RRT_object.collision_check_config(O,I, q_check)

%% plotting
figure(1) ; clf ; axis equal ; hold on ; grid on

plot(A)
plot(W)

if dimension == 3
    view(3)
end

%% returning a path

clear ; clc ;
goal = 5;
RRT_object = robot_arm_RRT_HLP_2D('goal', goal);

RRT_object.nodes = [1 2 3 4 5];
RRT_object.node_parent_idx = [3 0 4 2 1];
RRT_object.num_nodes = 5;

% path should return [2 4 3 1 5]

RRT_object.generate_waypoints()

path = RRT_object.waypoints;
check_path = [2 4 3 1 5]'
assert(isequal(path, check_path));

