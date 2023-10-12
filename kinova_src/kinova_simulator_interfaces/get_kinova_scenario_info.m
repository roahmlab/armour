function [scenario_obstacles,start,goal_radius,goal_type,goal] = get_kinova_scenario_info(scenario)

%%% WORLD PARAMETERS %%%
scenario_obstacles = {};
switch scenario
    case 1 % table
        start = [0; 0.5; 0; -0.5; 0; 0; 0];
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'
        goal = [0;-0.5;0; 0.5;0;0; 0];
        obs_center = [1.1; 0; 0.8];
        obs_side_lengths = [1, 4, 0.01];
%         obs_center = [0; 0; 1.3];
%         obs_side_lengths = [0.01, 4, 1];
        % create obstacle
        O = box_obstacle_zonotope('center',obs_center(:),...
                'side_lengths',obs_side_lengths) ;
            
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles({O});
        
    case 2 % wall/doorway
        start = [pi/2; 0.5; 0; 0; 0; 0;0];
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'
        goal = [-pi/2;0.5;0; 0.5;0;0;0];
        obs_center = [1.1; 0; 0.8];
        obs_side_lengths = [1, 0.01, 4];
%         obs_center = [0; 0; 1.3];
%         obs_side_lengths = [4, 0.01, 1];
        % create obstacle
        O = box_obstacle_zonotope('center',obs_center(:),...
                'side_lengths',obs_side_lengths) ;
            
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles({O});
    case 3 % posts
        start = [pi/2; pi/4; 0; 0; 0; 0; 0];
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'
        goal = [0.15;-0.75;0.2; 0.4;0.3;0.2; 0];
        obs_center = [0.8; -0.25; 2];
        obs_side_lengths = [0.05, 0.05, 4];
        % create obstacle
        O1 = box_obstacle_zonotope('center',obs_center(:),...
                'side_lengths',obs_side_lengths) ;
        obs_center = [0.4; 0.25; 2];
        obs_side_lengths = [0.05, 0.05, 4];
        % create obstacle
        O2 = box_obstacle_zonotope('center',obs_center(:),...
                'side_lengths',obs_side_lengths) ;
            
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles({O1, O2});        
        
    case 4 % shelves
        % manually create start
        start = [0;-0.5;0;0.5;0;0; 0] ; % on top of shelf 1 (in front of robot)

        % manually create goal
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'

        % configuration goals (set goal_type to 'configuration')
        % goal = [0;+1;0;-1;0;0] ; % bottom of front shelf
        % goal = [0.5;-0.5;0;0.5;0;0] ; % left of front shelf
        % goal = [0.5;+1;0;-1;0;0] ; bottom left of front shelf
        goal = [-pi/2 ; pi/2;-pi/2;0.5;0;0; 0] ; % top of left shelf
        % goal = [pi/2 ;+1;0;-1;0;0] ; % bottom of left shelf

        % end effector location goals (set goal_type to 'end_effector_location')
        % goal = [0; 1 ; 0.4] ; bottom middle of left shelf
        % goal = [0.25 ; 1 ; 1] ; top right of left shelf
        % goal = [1 ; -0.4 ; 0.6] ; % bottom right of front shelf
        % goal = [-0.3 ; 1 ; 0.6] ; % bottom of left shelf

        % shelf parameters
        shelf_center_1 = [1.1 ; 0 ; 0.7] ;
        shelf_center_2 = [0 ; 1.1 ; 0.7] ;
        shelf_height = 1.4 ; % m
        shelf_width = 1.2 ; % m 
        shelf_depth = 0.8 ; % m
        N_shelves = 3 ;
        min_shelf_height = 0.3 ;
        max_shelf_height = 1.3 ;
        
        % create shelves
        shelf_1 = make_shelf_obstacle(shelf_center_1,shelf_height,shelf_width,...
            shelf_depth,N_shelves,min_shelf_height,max_shelf_height,1) ;
        
        shelf_2 = make_shelf_obstacle(shelf_center_2,shelf_height,shelf_width,...
            shelf_depth,N_shelves,min_shelf_height,max_shelf_height,2) ;
        
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles([shelf_1, shelf_2]);
        
    case 5 % inside box
%         start = [0; -pi/12; 0; pi/3; 0; pi/3];
        start = [0; 0; 0; -pi/2; 0; 0; 0];
%         start = [0.2; 0.2; 0.2; 0.2; 0.2; 0.2];
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'
        goal = [0.15;0.1;0.2; 0.4;0.3;0.2; 0];
%         goal = [0.5; 0.2; -0.3; 0.4; -0.2; 0.7];
        
        total_box_side_lengths = [0.4, 0.4, 0.66];
        total_box_center = [0.45, 0, total_box_side_lengths(3)/2];
%         total_box_center = [0.7, 0, total_box_side_lengths(3)/2];

        
        obs_center = [total_box_center(1); total_box_center(2) + total_box_side_lengths(2)/2; total_box_center(3)];
        obs_side_lengths = [total_box_side_lengths(1), 0.01, total_box_side_lengths(3)];
        O1 = box_obstacle_zonotope('center',obs_center(:),'side_lengths',obs_side_lengths) ;
        
        obs_center = [total_box_center(1) - total_box_side_lengths(1)/2; total_box_center(2); total_box_center(3)];
        obs_side_lengths = [0.01, total_box_side_lengths(2), total_box_side_lengths(3)];
        O2 = box_obstacle_zonotope('center',obs_center(:),'side_lengths',obs_side_lengths) ;
        
        obs_center = [total_box_center(1); total_box_center(2) - total_box_side_lengths(2)/2; total_box_center(3)];
        obs_side_lengths = [total_box_side_lengths(1), 0.01, total_box_side_lengths(3)];
        O3 = box_obstacle_zonotope('center',obs_center(:),'side_lengths',obs_side_lengths) ;
        
        obs_center = [total_box_center(1) + total_box_side_lengths(1)/2; total_box_center(2); total_box_center(3)];
        obs_side_lengths = [0.01, total_box_side_lengths(2), total_box_side_lengths(3)];
        O4 = box_obstacle_zonotope('center',obs_center(:),'side_lengths',obs_side_lengths) ;
        
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles({O1, O2, O3, O4});
    case 6 % sink to cupboard
        start = [0; pi/6; 0; -pi/3-0.15; 0; -pi/3; 0];
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'
        goal = [pi/6;5*pi/12;-pi/2; -pi/8;pi/2;-pi/2; 0];
%       goal = [-pi/3;-pi/3;pi/6; pi/6;+pi/3;pi/6];

        % counter size
        counter_center = [0.6; 0; 0.6];
        counter_length = 0.5;
        counter_width = 2;
        
        sink_width = 0.5;
        sink_depth = 0.3;
        
        cupboard_center = [0.6; -0.55; 1.4];
        cupboard_length = counter_length;
        cupboard_width = 0.5;
        cupboard_depth = 0.5;
        
        % create counter
        obs_center = counter_center + [0; sink_width/2 + counter_width/2; 0];
        obs_side_lengths = [counter_length, counter_width, 0.01];
        O1 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = counter_center + [0; -sink_width/2 - counter_width/2; 0];
        obs_side_lengths = [counter_length, counter_width, 0.01];
        O2 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
                
        % create sink
        obs_center = counter_center + [0; sink_width/2; -sink_depth/2];
        obs_side_lengths = [sink_width, 0.01, sink_depth];
        O3 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = counter_center + [0; -sink_width/2; -sink_depth/2];
        obs_side_lengths = [sink_width, 0.01, sink_depth];
        O4 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = counter_center + [sink_width/2; 0; -sink_depth/2];
        obs_side_lengths = [0.01, sink_width, sink_depth];
        O5 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = counter_center + [-sink_width/2; 0; -sink_depth/2];
        obs_side_lengths = [0.01, sink_width, sink_depth];
        O6 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = counter_center + [0; 0; -sink_depth];
        obs_side_lengths = [sink_width, sink_width, 0.01];
        O7 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        % create cupboard
        obs_center = cupboard_center + [0; cupboard_width/2; 0];
        obs_side_lengths = [cupboard_length, 0.01, cupboard_depth];
        O8 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = cupboard_center + [0; -cupboard_width/2; 0];
        obs_side_lengths = [cupboard_length, 0.01, cupboard_depth];
        O9 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = cupboard_center + [0; 0; cupboard_depth/2];
        obs_side_lengths = [cupboard_length, cupboard_width, 0.01];
        O10 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = cupboard_center + [0; 0; -cupboard_depth/2];
        obs_side_lengths = [cupboard_length, cupboard_width, 0.01];
        O11 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        obs_center = cupboard_center + [cupboard_length/2; 0; 0];
        obs_side_lengths = [0.01, cupboard_width, cupboard_depth];
        O12 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles({O1, O2, O3, O4, O5, O6, O7, O8, O9, O10, O11, O12});
        
        
    case 7 % reach through window
%         start = [-pi/2; 0; pi/2; pi/4; 0; 0; 0];
        start = [0; pi/2; 0; -pi/4; 0; 0; 0];
        goal_radius = 0.05 ;
        goal_type = 'configuration' ; % 'configuration' or 'end_effector_location'
        goal = [0;0;0; 0;pi/3;pi/3; 0];
        
        % window size
        window_center = [0.6; 0; 0.8];
        window_side_length = 0.5;
        
        obs_height = 1.5;
        obs_width = 1.5;
        
        % create obstacle
        obs_center = window_center + [0; 0; -window_side_length/2 - obs_height/2];
        obs_side_lengths = [0.01, 4, obs_height];
        O1 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        % create obstacle
        obs_center = window_center + [0; 0; +window_side_length/2 + obs_height/2];
        obs_side_lengths = [0.01, 4, obs_height];
        O2 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        % create obstacle
        obs_center = window_center + [0; -window_side_length/2 - obs_width/2; 0];
        obs_side_lengths = [0.01, obs_width, 4];
        O3 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        % create obstacle
        obs_center = window_center + [0; +window_side_length/2 + obs_width/2; 0];
        obs_side_lengths = [0.01, obs_width, 4];
        O4 = box_obstacle_zonotope('center',obs_center(:),...
            'side_lengths',obs_side_lengths) ;
        
        scenario_obstacles = fetch_obstacles_to_kinova_obstacles({O1, O2, O3, O4});
        
        
    otherwise error('scenario not recognized');
end
end

function [o2] = fetch_obstacles_to_kinova_obstacles(o1)
    for i = 1:length(o1)
        obs_center = [o1{i}.center(3) - 0.8; o1{i}.center(2); o1{i}.center(1) + 0.25];
        side_lengths = [o1{i}.side_lengths(3), o1{i}.side_lengths(2), o1{i}.side_lengths(1)];
        o2{i} = box_obstacle_zonotope('center', obs_center(:), ...
            'side_lengths', side_lengths);
    end
end
%%% END WORLD PARAMETERS %%%
