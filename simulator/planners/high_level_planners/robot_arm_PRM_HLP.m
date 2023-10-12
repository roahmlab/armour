classdef robot_arm_PRM_HLP < robot_arm_sampling_based_HLP
    %% properties
    properties
        new_node_growth_distance = 0.1 ;
        new_node_nearest_neighbor_distance = 0.5 ;
    end
    
    %% methods
    methods
        %% constructor
        function HLP = robot_arm_PRM_HLP(varargin)
            HLP@robot_arm_sampling_based_HLP(varargin{:}) ;
        end
        
        %% setup
        function setup(HLP,agent_info,world_info)
            % call parent method
            setup@robot_arm_sampling_based_HLP(HLP,agent_info,world_info)
            
            % make sure the nearest neighbor distance is greater than the
            % new node growth distance
            if HLP.new_node_nearest_neighbor_distance < HLP.new_node_growth_distance
               HLP.vdisp('Updating PRM nearest neighbor distance!',5)
               HLP.new_node_nearest_neighbor_distance = 2.*HLP.new_node_growth_distance ;
            end
        end
        
        %% sampling
        function q_new = create_new_node(HLP,~,~)
            HLP.vdisp('Creating new node',9)
%             
%             % create a new node in the valid workspace of the arm
%             B = HLP.arm_joint_state_limits ;
%             q_new = rand_range(B(1,:),B(2,:))' ;
            
            % pick an existing node
            n = HLP.n_nodes ;
            q_cur_idx = round(rand_range(1,n)) ;
            q_cur = HLP.nodes(:,q_cur_idx) ;
            
            % create a new node by generating a random unit vector,
            % then growing in that direction from the current node
            u_rand = make_random_unit_vector(HLP.arm_n_links_and_joints) ;
            q_new = q_cur + u_rand.*HLP.new_node_growth_distance ;
            
            % make sure the new node lies within the arm's bounds
            q_new = bound_array_elementwise(q_new,...
                HLP.arm_joint_state_limits(1,:)',...
                HLP.arm_joint_state_limits(2,:)') ;
        end
        
        function add_node(HLP,q_new,~,~)
            q_dists = dist_point_to_points(q_new,HLP.nodes) ;
            q_near_log = q_dists <= HLP.new_node_nearest_neighbor_distance ;
            q_near_idxs = HLP.all_node_idxs(q_near_log) ;
            HLP.update_nodes_and_adjacency_matrix(q_new,q_near_idxs,...
                q_dists(q_near_log)) ;
        end
    end
end