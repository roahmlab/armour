classdef arm_end_effector_RRT_HLP < RRT_HLP
    properties
        start = zeros(3,1) ;
        buffer = 0 ;
        edge_feasibility_check_discretization = 0.01 ; 
        O_cell ; % cell array for obstacles
        current_waypoint_patch_data
    end
    
    methods
        %% constructor
        function HLP = arm_end_effector_RRT_HLP(varargin)
            HLP@RRT_HLP(varargin{:}) ;
            
            % set additional plot data
            HLP.plot_data.waypoint_arm_volume = [] ;
        end
        
        function setup(HLP,agent_info,world_info)
            % call superclass
            setup@RRT_HLP(HLP,agent_info,world_info) ;
            
            % set the start and goal
            J_start = agent_info.get_joint_locations(world_info.start) ;
            HLP.start = J_start(:,end) ;
            
            J_goal = agent_info.get_joint_locations(world_info.goal) ;
            HLP.goal = J_goal(:,end) ;
            
            % initialize tree
            HLP.initialize_tree(agent_info)
            
            % set the bounds from the buffer
            b = HLP.buffer ;
            HLP.bounds = world_info.bounds - b.*[-1 1 -1 1 -1 1] ; % 3D
        end
        
        function z = get_agent_position(HLP,agent_info)
            % get initial pose
            q = agent_info.state(agent_info.joint_state_indices,end) ;
            
            % get position
            HLP.vdisp('Getting starting end effector location',8)
            z_all = agent_info.get_joint_locations(q) ;
            z = z_all(:,end) ;
        end
        
        %% get waypoint
        function q = get_waypoint(HLP,agent_info,world_info,lookahead_distance)
            % get the waypoint in 3-D space
            z = get_waypoint@RRT_HLP(HLP,agent_info,world_info,lookahead_distance) ;
            HLP.current_waypoint = z ;
            
            % get the agent's current state
            q_0 = agent_info.state(agent_info.joint_state_indices,end) ;
            
            % inverse kinematics to get q close to z
            HLP.vdisp('Running inverse kinematics to get waypoint')
            [q,exitflag] = agent_info.inverse_kinematics(z,q_0) ;
            
            if exitflag < 0
                HLP.vdisp('IK failed! Using global goal config as waypoint',8)
                q = world_info.goal ;
            end
            
            % get patch data to plop
            HLP.current_waypoint_patch_data = agent_info.get_collision_check_volume(q) ;
        end
        
        %% tree growth
        function rand_node = sample(HLP,agent_info)
            % generate a random node within the bounds of the workspace
            if rand < HLP.goal_as_new_node_rate
                rand_node = HLP.goal ;
            else
                B = HLP.bounds ;
                switch HLP.dimension
                    case 2
                        rand_node = rand_range(B([1 3]),B([2 4]))' ;
                    case 3
                        rand_node = rand_range(B([1 3 5]),B([2 4 6]))' ;
                    otherwise
                        error('The HLP''s dimension must be 2 or 3!')
                end
            end
        end
        
        function exit_flag = grow_tree(HLP,agent_info,world_info)
            % update the HLP obstacles (this speeds up the edge feas check)
            O = world_info.obstacles ;
            b = HLP.buffer ;
            
            % get just the non-dynamic obstacles
            O_static = {} ;
            for idx = 1:length(O)
                o = O{idx} ;
                if isa(o,'box_obstacle_zonotope')
                    O_static = [O_static, {o}] ;
                end
            end
            
            % get the representation of the obstacles needed for the HLP
            if length(O_static) > length(HLP.O_cell)
                HLP.O_cell = cellfun(@(o) [o.center ; (2*b + o.side_lengths(:))],...
                    O_static,'UniformOutput',false) ;
            end
            
            % call the superclass method
            exit_flag = grow_tree@RRT_HLP(HLP,agent_info,world_info) ;
        end
        
        %% node feasibility check
        function out = edge_feasibility_check(HLP,node_A,node_B,~)
            % this function should return TRUE if the node is feasible
            
            % get the line between the two nodes
            l = [node_A, node_B] ;
            d_min = HLP.edge_feasibility_check_discretization ; % 1 cm
            d = vecnorm(l(:,2) - l(:,1)) ;
            
            if d > 0
                N = ceil(d/d_min) ;
                d_des = linspace(0,d,N) ;
                L = match_trajectories(d_des,[0 d],l) ;

                % check the obstacles
                O_chk = cellfun(@(o) dist_point_to_box(L,o(4),o(5),o(6),o(1:3)) == 0,...
                    HLP.O_cell,'UniformOutput',false) ;

                % check all the obstacles
                out = ~any(cell2mat(O_chk)) ;
            else
                out = false ;
            end
        end
        
        %% plotting
        function plot(HLP)
            % call superclass method
            plot@RRT_HLP(HLP)
            
            % plot the waypoint
            if HLP.plot_waypoint_flag && ~isempty(HLP.current_waypoint)
                if check_if_plot_is_available(HLP,'waypoint_arm_volume')
                    HLP.plot_data.waypoint_arm_volume.Faces = HLP.current_waypoint_patch_data.faces ;
                    HLP.plot_data.waypoint_arm_volume.Vertices = HLP.current_waypoint_patch_data.vertices ;
                else
                    HLP.plot_data.waypoint_arm_volume = patch(HLP.current_waypoint_patch_data,...
                        'FaceColor','g','FaceAlpha',0.1) ;
                end
            end
        end
    end
end