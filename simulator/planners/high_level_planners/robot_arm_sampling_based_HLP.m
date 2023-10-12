classdef robot_arm_sampling_based_HLP < high_level_planner
    %% properties
    properties
        % arm
        arm_n_states
        arm_n_inputs
        arm_n_links_and_joints
        arm_joint_state_limits
        arm_joint_speed_limits
        arm_joint_input_limits
        arm_joint_state_indices
        arm_joint_speed_indices
        
        % samples
        nodes
        adjacency_matrix
        n_nodes
        n_nodes_max = 20000 ;
        all_node_idxs ;
        sampling_timeout = 0.1 ; % seconds per planning iteration
        choose_goal_as_new_node_frequency = 0.1 ;
        make_new_graph_every_iteration_flag = true ;
        
        % path
        start
        best_path_nodes
        best_path_node_idxs
        
        % plotting
        plot_while_sampling_flag = false ;
        plot_HLP_flag = false ;
        plot_waypoint_color = [0 1 0] ;
    end
    
    %% methods
    methods
        %% constructor
        function HLP = robot_arm_sampling_based_HLP(varargin)
            % robot_arm_sampling_based_HLP < high_level_planner
            %
            % This is the superclass for sampling-based high-level planners
            % to be used with the robot arms in the simulator framework.
            
            % create the planner
            HLP@high_level_planner(varargin{:}) ;
            HLP.all_node_idxs = 1:HLP.n_nodes_max ;
            
            % set up the plot data
            HLP.plot_data.new_node = [] ;
            HLP.plot_data.waypoint = [] ;
        end
        
        function setup(HLP,agent_info,world_info)
            % get all the necessary arm properties filled in
            HLP.vdisp('Filling in HLP''s arm properties',9)
            HLP = fill_in_arm_properties(HLP,agent_info,true) ;
            
            % set the HLP's dimension
            HLP.dimension = agent_info.dimension ;
            
            % get the world goal
            HLP.vdisp('Setting goal',9)
            HLP.start = world_info.start ;
            HLP.goal = world_info.goal ;
            
            % set up graph structure
            HLP.initialize_graph(agent_info,world_info)
        end
        
        function initialize_graph(HLP,agent_info,world_info)
            % initialize nodes and adjacency matrix
            HLP.vdisp('Initializing planner graph structures',4)
            HLP.n_nodes = 1 ;
            HLP.nodes = nan(HLP.arm_n_links_and_joints,HLP.n_nodes_max) ;
            HLP.adjacency_matrix = sparse(HLP.n_nodes_max, HLP.n_nodes_max) ;
            
            % start the graph at the agent's current state (if the graph is
            % re-initialized at the beginning of each planning iteration,
            % then this ensures the graph starts where the agent is --
            % i.e., when HLP.make_new_graph_every_iteration_flag is true)
            HLP.nodes(:,1) = agent_info.state(HLP.arm_joint_state_indices,end) ;
            
            % reset best path and best path nodes
            HLP.best_path_nodes = [] ;
            HLP.best_path_node_idxs = [] ;
        end
        
        %% get waypoint
        function waypoint = get_waypoint(HLP,agent_info,world_info,lookahead_distance)
            % waypoint = HLP.get_waypoint(agent_info,world_info,lookahead_distance)
            %
            % Given the agent and world info, and a desired lookahead
            % distance, sample to generate a graph, then generate a
            % waypoint along the best path found on the graph.
            
            % run sampling algorithm to grow graph
            HLP.sample(agent_info,world_info)
            
            % get best path in graph
            HLP.find_best_path(agent_info,world_info) ;
            
            % generate waypoint
            if nargin < 4
                lookahead_distance = HLP.default_lookahead_distance ;
            end
            waypoint = HLP.create_waypoint(agent_info,world_info,lookahead_distance) ;
            
            if HLP.plot_HLP_flag
                HLP.vdisp('Plotting waypoint',6)
                
                V_plot = agent_info.get_collision_check_volume(waypoint) ;
                
                if check_if_plot_is_available(HLP,'waypoint')
                    HLP.plot_data.waypoint.Faces = V_plot.faces ;
                    HLP.plot_data.waypoint.Vertices = V_plot.vertices ;
                else
                    HLP.plot_data.waypoint = patch(V_plot,...
                        'facecolor',HLP.plot_waypoint_color,'facealpha',0.1) ;
                end
                drawnow()
            end
        end
        
        function waypoint = create_waypoint(HLP,agent_info,world_info,lookahead_distance)
            % waypoint = HLP.create_waypoint(agent_info,world_info,lookahead_distance)
            %
            % Given the HLP's current best path, get the cumulative
            % distance along the path, then interpolate by the lookahead
            % distance to create a waypoint along the path.
            
            HLP.vdisp('Making waypoint',5) ;
            
            % get the best path
            path_nodes = HLP.best_path_nodes ;
            
            % if the best path isn't empty, get the distance along it
            if ~isempty(path_nodes)
                % get the current state
                z_cur = agent_info.state(:,end) ;
                q_cur = z_cur(HLP.arm_joint_state_indices) ;
                
                % get the distance of the current node along the path
                [d,~,d_along] = dist_point_on_polyline(q_cur,path_nodes) ;
                
                % look ahead to generate a waypoint
                d_lkhd = min([d + lookahead_distance, d_along(end)]) ;
                waypoint = match_trajectories(d_lkhd,d_along,path_nodes) ;
                
%                 % get the first waypoint that is further along than the
%                 % current state
%                 q_idx = find(d_along > d,1,'first') ;
%                 
%                 if ~isempty(q_idx)
%                     waypoint = path_nodes(:,q_idx) ;
%                 else
%                     waypoint = path_nodes(:,end) ;
%                 end
            else
                waypoint = HLP.goal ;
            end
        end
        
        %% graph search
        function find_best_path(HLP,agent_info,world_info)
            % HLP.find_best_path(agent_info,world_info)
            %
            % Given the current agent and world information, find the best
            % path to the goal. By default, this finds the shortest path
            % (in the graph represented by HLP.adjacency_matrix) from the
            % node closest to the agent's current state to the node closest
            % to the global goal position.
            
            HLP.vdisp('Finding best path',5)
            
            % find the node closest to the goal state
            [~,q_goal_idx] = HLP.find_nearest_node(HLP.goal) ;
            
            % find the best path between the start and goal nodes
            [~,best_path_idxs,~] = graphshortestpath(HLP.adjacency_matrix,...
                1,q_goal_idx) ;
            
            HLP.best_path_nodes = HLP.nodes(:,best_path_idxs) ;
            HLP.best_path_node_idxs = best_path_idxs ;
        end
        
        function [q_near,q_near_idx] = find_nearest_node(HLP,q)
            q_dists = dist_point_to_points(q,HLP.nodes) ;
            [~,q_near_idx] = min(q_dists,[],'omitnan') ;
            q_near = HLP.nodes(:,q_near_idx) ;
            
            if nargout < 2
                q_near_idx = [] ;
            end
        end
        
        %% sampling
        function sample(HLP,agent_info,world_info)
            % prep for while loop
            start_tic = tic ;
            t_cur = toc(start_tic) ;
            
            if HLP.make_new_graph_every_iteration_flag
                HLP.initialize_graph(agent_info,world_info)
            end
            
            % run sampling
            HLP.vdisp('Running sampling algorithm',5)
            
            while HLP.n_nodes <= HLP.n_nodes_max && t_cur <= HLP.sampling_timeout
                % generate a new sample
                q_new = HLP.create_new_node(agent_info,world_info) ;
                
                % check if the new node is feasible
                q_is_feasible = HLP.check_if_node_is_feasible(q_new,agent_info,world_info) ;
                
                % if the new node is not in collision, add it to graph
                if q_is_feasible
                    HLP.vdisp('Feasible node found!',9)
                    HLP.add_node(q_new,agent_info,world_info) ;
                    new_node_color = [0 0 1] ;
                else
                    HLP.vdisp('Node infeasible',9)
                    new_node_color = [1 0 0] ;
                end
                
                % plot the new node (FOR DEBUGGING)
                if HLP.plot_while_sampling_flag
                    V_plot = agent_info.get_collision_check_volume(q_new) ;
                    if check_if_plot_is_available(HLP,'new_node')
                        HLP.plot_data.new_node.Faces = V_plot.faces ;
                        HLP.plot_data.new_node.Vertices = V_plot.vertices ;
                    else
                        HLP.plot_data.new_node = patch(V_plot,'facecolor',new_node_color,'facealpha',0.05) ;
                    end
                    drawnow()
                end
                
                % update timing
                t_cur = toc(start_tic) ;
            end
            HLP.vdisp(['Number of nodes total: ',num2str(HLP.n_nodes)],5)
            HLP.vdisp(['Time spent: ',num2str(t_cur,'%0.2f'),' seconds'],8)
        end
        
        %% sample generation
        function q_new = create_new_node(HLP,agent_info,world_info)
            HLP.vdisp('Sampling a new node',8)
            q_new = rand_range(HLP.arm_joint_state_limits(1,:)',...
                HLP.arm_joint_state_limits(2,:)') ;
        end
        
        %% feasibility checking
        function out = check_if_node_is_feasible(HLP,q_to_check,agent_info,world_info)
            HLP.vdisp('Checking if node is feasible',9)
            out = ~(HLP.collision_check_single_state(q_to_check,agent_info,world_info)) ;
        end
        
        function out = collision_check_single_state(HLP,q_to_check,agent_info,world_info)
            HLP.vdisp('Collision checking node',9)
            O = world_info.obstacles ;
            N_O = length(O) ;
            out = false ;
            o_idx = 1 ;
            V_arm = agent_info.get_collision_check_volume(q_to_check) ;
            
            while (o_idx <= N_O) && ~out
                O_idx = O{o_idx} ;
                out = HLP.collision_check_single_obstacle(O_idx,V_arm) ;
                o_idx = o_idx + 1 ;
            end
        end
        
        function out = collision_check_single_obstacle(HLP,obstacle_object,arm_volume)
            HLP.vdisp('Collision checking obstacle',10)
            switch HLP.dimension
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
        
        %% node and graph updates
        function add_node(HLP,q_new,agent_info,world_info)
            % HLP.add_node(q_new,agent_info,world_info)
            %
            % Add the node q_new to the list of nodes, and update the graph
            % representation. By default, the new node is connected to its
            % nearest neighbor.
            
            [~,q_near_idx] = HLP.find_nearest_node(q_new) ;
            HLP.update_nodes_and_adjacency_matrix(q_new,q_near_idx) ;
        end
        
        function update_nodes_and_adjacency_matrix(HLP,q_new,q_connect_idxs,weights)
            % update_nodes_and_adjacency_matrix(HLP,q_new,q_connect_idxs)
            %
            % Increment the number of nodes, add the new node to the node
            % list, and update the adjacency matrix so that the new node is
            % connected to all nodes indexed by q_connect_idxs.
            
            % update the number of nodes
            HLP.n_nodes = HLP.n_nodes + 1 ;
            
            % add the node to the list of nodes
            HLP.nodes(:,HLP.n_nodes) = q_new ;
            
            % update the adjacency matrix
            if nargin < 4
                weights = 1 ;
            end
            HLP.adjacency_matrix(q_connect_idxs,HLP.n_nodes) = weights ;
            HLP.adjacency_matrix(HLP.n_nodes,q_connect_idxs) = weights ;
        end
        
        %% utility
        function out = check_arm_properties(HLP)
            % out = HLP.check_arm_properties()
            %
            % Checks each property that starts with 'arm' and makes sure
            % they are all filled in by the planner that is using this HLP.
            
            all_property_names = fieldnames(HLP) ;
            
            idx = 1 ;
            out = true ;
            
            while idx <= length(all_property_names) && out
                current_property = all_property_names{idx} ;
                
                if strcmpi(current_property(1:3),'arm') && isempty(HLP.(current_property))
                    out = false ;
                    HLP.vdisp(['Missing property: ',current_property],3) ;
                end
                
                idx = idx + 1 ;
            end
        end
    end
end