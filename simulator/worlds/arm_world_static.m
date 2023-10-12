classdef arm_world_static < world
    properties
        % setup info
        include_base_obstacle = true ;
        obstacle_size_range = [0.05 0.15] ; % [min, max] side length
        create_configuration_timeout = 1 ;
        create_obstacle_timeout =  1 ;
        min_dist_in_config_space_between_start_and_goal
        floor_normal_axis = 3;
        goal_type = 'configuration';
        
        % arm info
        arm_joint_state_limits
        arm_n_links_and_joints
        arm_joint_state_indices
        
        % goal plotting
        goal_plot_patch_data
        goal_plot_face_color = [0 1 0] ;
        goal_plot_face_alpha = 0.1 ;
        goal_plot_edge_color = [0 1 0] ;
        goal_plot_edge_alpha = 0.3 ;
        goal_plot_edge_style = '--' ;
    end
    
    methods
        %% constructor
        function W = arm_world_static(varargin)
            % W = arm_world_static('Aproperty1',value1,'property2',value2,...)
            
            default_goal_radius = 0.05 ; % rad/joint
            W@world('start',[],'goal',[],'N_obstacles',0,...
                'goal_radius',default_goal_radius,...
                varargin{:}) ;
            
            W.plot_data.obstacles = [] ;
            W.plot_data.goal = [] ;
        end
        
        %% setup
        function setup(W,I)
            W.vdisp('Running arm world setup',1)
            
            % make sure the world's properties agree with the arm's
            W.get_arm_info(I)
            
            % generate an obstacle to attach the arm's baselink to
            if W.include_base_obstacle
                W.obstacles = {W.create_base_obstacle()} ;
            end
            
            % create a start configuration
            if isempty(W.start)
                W.create_start(I) ;
            end
            
            % create obstacles
            if isempty(W.obstacles) || ...
               (W.include_base_obstacle && length(W.obstacles) == 1)
                for idx = 1:W.N_obstacles
                    O = W.create_collision_free_obstacle(I) ;
                    if ~isempty(O)
                        W.obstacles = [W.obstacles, {O}] ;
                    end
                end
            end
            
            % update the number of obstacles
            W.N_obstacles = length(W.obstacles) ;
            
            % create a goal configuration
            if isempty(W.goal)
                W.create_goal(I) ;
            end
            W.goal_plot_patch_data = I.get_collision_check_volume(W.goal) ;
            
            W.vdisp('Arm world setup complete',2)
        end
        
        %% get bounds and joint limits
        function get_arm_info(W,I)
            % update the world's dimension based on the arm's dimension
            if W.dimension ~= I.dimension
                W.vdisp(['Updating world dimension to ',num2str(I.dimension),...
                    ' to match the arm!'],3)
                W.dimension = I.dimension ;
            end
            
            % set world bounds based on agent limits
            W.bounds = I.reach_limits ;
            
            % set any joint limits that are +Inf to pi and -Inf to -pi
            joint_state_limits = I.joint_state_limits ;
            joint_limit_infs = isinf(joint_state_limits) ;
            joint_state_limits(1,joint_limit_infs(1,:)) = -pi ;
            joint_state_limits(2,joint_limit_infs(2,:)) = +pi ;
            
            W.arm_joint_state_limits = joint_state_limits ;
            W.arm_n_links_and_joints = size(joint_state_limits,2) ;
            W.arm_joint_state_indices = I.joint_state_indices ;
            
            % set minimum distance between start and goal based on the
            % joint limits
            joint_ranges = diff(joint_state_limits,[],1) ;
            W.min_dist_in_config_space_between_start_and_goal = norm(0.25*joint_ranges) ;
        end
        
        %% make start and goal
        function create_start(W,I)
            W.vdisp('Making start configuration',5)
            W.start = W.create_collision_free_configuration(I) ;
            if isempty(W.start)
                W.vdisp('Using agent current pose as start config',3)
                W.start = I.state(I.joint_state_indices,end) ;
            end
        end
        
        function create_goal(W,I)
            W.vdisp('Making goal configuration',5)
            
            dist_between_start_and_goal = 0 ;
            start_tic = tic ;
            t_cur = toc(start_tic) ;
            new_goal = [] ;
            
            while dist_between_start_and_goal < W.min_dist_in_config_space_between_start_and_goal && ...
                t_cur <= W.create_configuration_timeout
            
                % new_goal = rand_range(W.arm_joint_state_limits(1,:),W.arm_joint_state_limits(2,:))' ;
                new_goal = W.create_collision_free_configuration(I) ;
                
                dist_between_start_and_goal = norm(W.start - new_goal) ;
                
                t_cur = toc(start_tic) ;
            end
            
            if isempty(new_goal)
                W.vdisp('Goal creation failed! Using random goal',3)
                W.goal = rand_range(W.arm_joint_state_limits(1,:),W.arm_joint_state_limits(2,:))' ;
            else
                W.goal = new_goal ;
            end
            
        end
        
        %% make configurations
        function q = create_collision_free_configuration(W,I)
            config_is_valid = false ;
            start_tic = tic ;
            t_cur = toc(start_tic) ;
            
            while ~config_is_valid && t_cur <= W.create_configuration_timeout
                q = W.create_random_configuration() ;
                config_is_valid = ~(W.collision_check_single_state(I,q)) ;
                t_cur = toc(start_tic) ;
            end
            
            if ~config_is_valid
                q = [] ;
                W.vdisp('Configuration creation failed!',3)
            end
        end
        
        function q = create_random_configuration(W)
            q = rand_range(W.arm_joint_state_limits(1,:),W.arm_joint_state_limits(2,:))' ;
        end
        
        %% make obstacles
        function O = create_collision_free_obstacle(W,I,q)
            if nargin < 3
                q = W.start ;
            end
            
            obstacle_is_valid = false ;
            start_tic = tic ;
            t_cur = toc(start_tic) ;
            V_arm = I.get_collision_check_volume(q) ;
            
            while ~obstacle_is_valid && t_cur <= W.create_obstacle_timeout
                O = W.create_random_obstacle() ;
                obstacle_is_valid = ~(W.collision_check_single_obstacle(O,V_arm)) ;
                t_cur = toc(start_tic) ;
            end
            
            if ~obstacle_is_valid
                O = [] ;
                W.vdisp('Obstacle creation failed! Try again...',3)
            end
        end
        
        function O = create_random_obstacle(W)
            % create center
            B = W.bounds ;
            center = [rand_range(B(1),B(2)) ; rand_range(B(3),B(4))] ;
%             center = [rand_range(0,B(2)) ; rand_range(B(3),B(4))] ;
            
            if W.dimension == 3
                center = [center ; rand_range(B(5),B(6))] ;
            end
            
            % create side lengths
            side_lengths = rand_range(W.obstacle_size_range(1),...
                W.obstacle_size_range(2),[],[],1,W.dimension) ;
            
            % create obstacle
            O = box_obstacle_zonotope('center',center(:),...
                'side_lengths',side_lengths) ;
        end
        
        
        function O = create_base_obstacle(W)
            W.vdisp('Making base obstacle',3) ;
            switch W.dimension
                case 2
                    side_lengths = [W.bounds(2) - W.bounds(1), 0.05] ;
                    center = [0 ; -side_lengths(2)/2] ;
                case 3
                    switch W.floor_normal_axis
                        case 1
                            side_lengths = [0.05, ...
                                W.bounds(4) - W.bounds(3), ...
                                W.bounds(6) - W.bounds(5)] ;
                            center = [-side_lengths(1)/2 ; 0 ; 0] ;
                        case 2
                            side_lengths = [W.bounds(2) - W.bounds(1), ...
                                0.05, ...
                                W.bounds(6) - W.bounds(5)] ;
                            center = [0 ; -side_lengths(2)/2 ; 0] ;
                        case 3
                            side_lengths = [W.bounds(2) - W.bounds(1), ...
                                W.bounds(4) - W.bounds(3),...
                                0.05] ;
                            center = [0 ; 0 ; -side_lengths(3)/2] ;
                        otherwise
                            error('floor axis not supported');
                    end
            end
            
            O = box_obstacle_zonotope('center',center,...
                'side_lengths',side_lengths,...
                'plot_face_color',[0.1 0.1 0.5],...
                'plot_edge_color',[0 0 1]) ;
        end
        
        %% collision checking
        function out = collision_check_single_state(W,I,q)
            % out = collision_check_single_state(W,agent_info,agent_state)
            %
            % Run a collision check for the given state and return true if
            % it is in collision. This gets called by W.collision_check.
            
            O = W.obstacles ;
            N_O = length(O) ; % in case W.N_obstacles is wrong
            out = false ; % optimism!
            o_idx = 1 ;
            V_arm = I.get_collision_check_volume(q) ;
            
            while (o_idx <= N_O) && ~out
                O_idx = O{o_idx} ;
                out = W.collision_check_single_obstacle(O_idx,V_arm) ;
                o_idx = o_idx + 1 ;
            end
        end
        
        function out = collision_check_single_obstacle(W,obstacle_object,arm_volume)
            switch W.dimension
                case 2
                    obstacle_object = obstacle_object.collision_check_patch_data.vertices ;
                    obstacle_object = [obstacle_object ; obstacle_object(1,:)] ;
                    [x_int,~] = polyxpoly(arm_volume(1,:)',arm_volume(2,:)',...
                        obstacle_object(:,1),obstacle_object(:,2)) ;
                    if isempty(x_int)
                        out = false ;
                    else
                        out = true ;
                    end
                case 3
                    O_str.faces = obstacle_object.collision_check_patch_data.faces ;
                    O_str.vertices = obstacle_object.collision_check_patch_data.vertices ;
                    check = SurfaceIntersection(O_str,arm_volume) ;
                    out = any(check(:)) ;
            end
        end
        
        %% goal check
        function out = goal_check(W,agent_info)
            
            z = agent_info.state(W.arm_joint_state_indices,:) ;
            
            switch W.goal_type
                case 'configuration'
                    dz = abs(z - repmat(W.goal,1,size(z,2))) ;
                    dz_log = dz <= W.goal_radius ;
                    out = any(all(dz_log,1)) ;
                case 'end_effector_location'
                    % get the joint locations
                    J = agent_info.get_joint_locations(z) ;
                    
                    if ~iscell(J)
                        J = {J};
                    end
                    
                    % concatenate all the end effector locations into one
                    % array
                    N_J = length(J) ;
                    J_ee = nan(3,N_J) ;
                    for idx = 1:N_J
                        J_ee(:,idx) = J{idx}(:,end) ;
                    end
                    
                    % check how far the end effector is from the goal
                    % location
                    dz = vecnorm(J_ee - repmat(W.goal_in_workspace,1,N_J)) ;
                    out = any(dz <= W.goal_radius) ;
                otherwise
                    error(['The goal type ',W.goal_type,' is not supported!'])
            end
            
        end
        
        %% plotting
        function plot(W)
            hold_check = ~ishold ;
            if hold_check
                hold on ;
            end
            
            % plot obstacles (each obstacle takes care of its own plotting)
            W.vdisp('Plotting obstacles',6)
            for idx = 1:W.N_obstacles
                plot(W.obstacles{idx}) ;
            end
            
            % plot goal config
            W.plot_goal()
            
            if hold_check
                hold off
            end
        end
        
        function plot_goal(W)
            W.vdisp('Plotting goal',6)
            
            hold_check = ~ishold ;
            if hold_check
                hold on ;
            end
            
            G = W.goal_plot_patch_data ;
            
            if ~check_if_plot_is_available(W,'goal')
                switch W.dimension
                    case 2
                        data = plot(G(1,:),G(2,:),'Color',W.goal_plot_edge_color,...
                            'LineStyle',W.goal_plot_edge_style) ;
                    case 3
                        data = patch(G,...
                            'LineStyle',W.goal_plot_edge_style,...
                            'FaceColor',W.goal_plot_face_color,...
                            'FaceAlpha',W.goal_plot_face_alpha,...
                            'EdgeColor',W.goal_plot_edge_color,...
                            'EdgeAlpha',W.goal_plot_edge_alpha) ;
                end
                
                W.plot_data.goal = data ;
            end
            
            if hold_check
                hold off
            end
        end
    end
end