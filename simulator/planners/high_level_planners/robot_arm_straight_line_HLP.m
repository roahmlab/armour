classdef robot_arm_straight_line_HLP < high_level_planner
    %% properties
    properties
        make_new_graph_every_iteration_flag = false ;
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
        
        % path
        best_path_nodes
        best_path_node_idxs
    end
    methods
        %%  constructor
        function HLP = robot_arm_straight_line_HLP(varargin)
            HLP@high_level_planner(varargin{:}) ;
        end
        
        %% setup
        function setup(HLP,agent_info,world_info)
            % get all the necessary arm properties filled in
            HLP.vdisp('Filling in HLP''s arm properties',9)
            HLP = fill_in_arm_properties(HLP,agent_info,false) ;
            
            % get the world goal
            HLP.vdisp('Setting goal',9)
            HLP.goal = world_info.goal ;
        end
        
        %% get waypoint
        function waypoint = get_waypoint(HLP,agent_info,~,lookahead_distance)
            q_cur = agent_info.state(HLP.arm_joint_state_indices, end) ;
            q_goal = HLP.goal ;
            dir_des = q_goal - q_cur ;
            full_rotate_joints = agent_info.joint_state_limits(1,:) == -Inf;
            dir_des(full_rotate_joints) = angdiff(q_cur(full_rotate_joints), ...
                                                  q_goal(full_rotate_joints));
            dir_des = dir_des./norm(dir_des) ;
            if nargin < 4
                lookahead_distance = HLP.default_lookahead_distance ;
            end
            waypoint = q_cur + lookahead_distance.*dir_des ;
        end
    end
end