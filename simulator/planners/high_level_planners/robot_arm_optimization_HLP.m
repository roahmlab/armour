classdef robot_arm_optimization_HLP < high_level_planner
    %% properties
    properties
        % arm
        arm_dimension
        arm_n_states
        arm_n_inputs
        arm_n_links_and_joints
        arm_joint_axes
        arm_joint_locations
        arm_joint_state_limits
        arm_joint_speed_limits
        arm_joint_input_limits
        arm_joint_state_indices
        arm_joint_speed_indices
        
        % optimization params
        buffer_dist = 0.1 ; % m
        get_new_waypoint_every_iteration_flag = false ;
        
        % plotting
        waypoint_plot_input_data
    end
    
    %% methods
    methods
        %% constructor
        function HLP = robot_arm_optimization_HLP(varargin)
            % create the planner
            HLP@high_level_planner(varargin{:}) ;
            
            % set up the plot data
            HLP.plot_data.waypoint = [] ;
        end
        
        function setup(HLP,agent_info,world_info)
            % get all the necessary arm properties filled in
            HLP.vdisp('Filling in HLP''s arm properties',9)
            HLP = fill_in_arm_properties(HLP,agent_info,true) ;
            
            % set the HLP's dimension
            HLP.dimension = agent_info.dimension ;
            
            % get the world goal
            HLP.vdisp('Getting world goal',9)
            HLP.goal = world_info.goal ;
            
            % set the initial waypoint
            if ~HLP.get_new_waypoint_every_iteration_flag
                HLP.vdisp('Running optimization!',6)
                
                % get start and goal
                x_start = world_info.start ;
                x_goal = HLP.goal ;
                
                % try running optimization
                [x_opt,exitflag] = HLP.run_optimization(x_start,x_goal,world_info) ;
                
                % if successful, set the new waypoint as the current
                % waypoint
                if exitflag > 0
                    HLP.vdisp('New waypoint found!',4)
                    HLP.current_waypoint = x_opt ;
                    HLP.waypoint_plot_input_data = agent_info.get_collision_check_volume(x_opt) ;
                end
            else
                exitflag = -1 ;
            end
            
            if exitflag <= 0
                HLP.vdisp('Setting initial waypoint as global goal',9)
                HLP.current_waypoint = HLP.goal ;
                HLP.waypoint_plot_input_data = agent_info.get_collision_check_volume(HLP.goal) ;
            end
        end
        
        %% get waypoint
        function waypoint = get_waypoint(HLP,agent_info,world_info,~)
            if HLP.get_new_waypoint_every_iteration_flag
                HLP.vdisp('Running optimization!',6)
                
                % get start and goal
                x_start = agent_info.state(HLP.arm_joint_state_indices,end) ;
                x_goal = HLP.goal ;
                
                % try running optimization
                [x_opt,exitflag] = HLP.run_optimization(x_start,x_goal,world_info) ;
                
                % if successful, set the new waypoint as the current
                % waypoint
                if exitflag > 0
                    HLP.vdisp('New waypoint found!',4)
                    HLP.current_waypoint = x_opt ;
                    HLP.waypoint_plot_input_data = agent_info.get_collision_check_volume(x_opt) ;
                end
            end
            
            waypoint = HLP.current_waypoint ;
        end
        
        %% run optimization
        function [x_opt,exitflag] = run_optimization(HLP,x_start,x_goal,world_info)
            obs_bounds = HLP.process_obstacles(world_info) ;
            
            % get start and goal end effector locations
            J_start = HLP.get_collision_check_locations(x_start) ;
            J_start = J_start(:,end) ;
            J_goal = HLP.get_collision_check_locations(x_goal) ;
            J_goal = J_goal(:,end) ;
            
            % make cost and constraints
            cost = @(x) HLP.cost_function(x,J_start,J_goal) ;
            nonlcon = @(x) HLP.nonlcon_function(x,obs_bounds) ;
            
            % make initial guess
            x_0 = 0.5.*(x_start + x_goal) ;
            
            % call fmincon
            [x_opt,~,exitflag] = fmincon(cost,x_0,[],[],[],[],...
                HLP.arm_joint_state_limits(1,:)',...
                HLP.arm_joint_state_limits(2,:)',...
                nonlcon) ;
        end
        
        %% optimization helper functions
        function c = cost_function(HLP,x,J_start,J_goal)
            J = HLP.get_collision_check_locations(x) ;
            c = vecnorm(J_goal - J(:,end)).^2 + vecnorm(J_start - J(:,end)).^2 ;
        end
        
        function [n,neq] = nonlcon_function(HLP,x,obs_bounds)
            J = get_collision_check_locations(HLP,x) ;
            n = -(cell2mat(cellfun(@(B) dist_point_to_box(J,B),...
                obs_bounds,'UniformOutput',false)) - HLP.buffer_dist);
            n = n(:) ;
            neq = [] ;
        end
        
        function obs_bounds = process_obstacles(HLP,world_info)
            obs_bounds = [] ;
            for idx = 1:length(world_info.obstacles)
                o = world_info.obstacles{idx} ;
                if isa(o,'box_obstacle_zonotope')
                    obs_bounds = [obs_bounds, {box_to_bounds(o.side_lengths(1),...
                        o.side_lengths(2),o.side_lengths(3),...
                        o.center)}] ;
                end
            end
        end
        
        function J = get_collision_check_locations(HLP,q)
            j_vals = q ; % joint angles
            j_locs = HLP.arm_joint_locations ; % joint locations
            
            % extract dimensions
            n = HLP.arm_n_links_and_joints ;
            d = HLP.arm_dimension ;
            
            % set up translations and rotations
            R_pred = eye(d) ;
            T_pred = zeros(d,1) ;
            
            % allocate array for the joint locations
            J = nan(d,n) ;
            
            % move through the kinematic chain and get the rotations and
            % translation of each link
            for idx = 1:n
                % get the value and location of the current joint
                j_idx = j_vals(idx) ;
                j_loc = j_locs(:,idx) ;
                
                % rotation matrix of current joint
                axis_pred = R_pred*HLP.arm_joint_axes(:,idx) ;
                R_succ = axis_angle_to_rotation_matrix_3D([axis_pred', j_idx])*R_pred ;
                
                % create translation
                T_succ = T_pred + R_pred*j_loc(1:d) - R_succ*j_loc(d+1:end) ;
                
                % fill in the joint location
                j_loc_local = j_locs((d+1):end,idx) ;
                J(:,idx) = -R_succ*j_loc_local + T_succ ;
                
                % update predecessors for next iteration
                R_pred = R_succ ;
                T_pred = T_succ ;
            end
            
            % add more points (hardcoded for now to be 3x the number of
            % joints)
            t = 1:n ;
            t2 = linspace(1,n,3*n) ;
            J = match_trajectories(t2,t,J) ;
        end
        
        %% plotting
        function plot(HLP)
            V = HLP.waypoint_plot_input_data ;
            
            if isempty(V)
                warning('Please run the planner/HLP setup method!')
            else
                if check_if_plot_is_available(HLP,'waypoint')
                    HLP.plot_data.waypoint.Faces = V.faces ;
                    HLP.plot_data.waypoint.Vertices = V.vertices ;
                else
                    HLP.plot_data.waypoint = patch(V,'facecolor',[0.5 1 0.5],'facealpha',0.2) ;
                end
            end
        end
    end
end