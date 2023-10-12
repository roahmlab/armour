classdef robot_arm_generic_planner < planner
    %% properties
    properties
        % arm
        dimension
        arm_n_states
        arm_n_inputs
        arm_n_links_and_joints
        arm_joint_state_limits
        arm_joint_speed_limits
        arm_joint_input_limits
        arm_joint_state_indices
        arm_joint_speed_indices
        
        % world
        start
        goal
        
        % planner
        t_stop = 1 ; % duration of stopping trajectory
        lookahead_distance = 1 ; % by 2-norm in config space
    end
    
    %% methods
    methods
        %% constructor
        function P = robot_arm_generic_planner(varargin)
            % create default high-level planner
            HLP = robot_arm_RRT_HLP() ;
            
            P@planner('HLP',HLP,varargin{:}) ;
            
            P.HLP.make_new_graph_every_iteration_flag = 1 ;
            P.HLP.sampling_timeout = P.t_plan ;
            
            % check that t_stop is greater than t_move
            if P.t_stop <= P.t_move
                t_stop_new = P.t_move + 1 ;
                P.vdisp(['Setting stopping time to ',num2str(t_stop_new),' seconds'],3)
                P.t_stop = t_stop_new ;
            end
        end
        
        %% setup
        function setup(P,agent_info,world_info)
            % fill in agent info
            P.vdisp('Getting agent info',9)
            P.dimension = agent_info.dimension ;
            P.bounds = agent_info.reach_limits ;
            P = fill_in_arm_properties(P,agent_info,false) ;
            
            % fill in world info
            P.vdisp('Getting world info',9)
            P.start = world_info.start ;
            P.goal = world_info.goal ;
            
            % set up high-level planner
            P.HLP.verbose = P.verbose ;
            P.HLP.setup(agent_info,world_info) ;
        end
        
        %% replan
        function [T_ref,U_ref,Z_ref] = replan(P,agent_info,world_info)
            % call HLP sample method
            P.HLP.sample(agent_info,world_info)
            
            % get HLP's best path
            P.HLP.find_best_path(agent_info,world_info) ;
            Z = P.HLP.best_path_nodes ;
            N_Z = size(Z,2) ;
            
            % generate time vector
            T_ref = linspace(0,P.t_move,N_Z) ;
            
            % generate joint speeds
            Zd = [agent_info.state(P.arm_joint_speed_indices,end),...
                diff(Z,[],2)./diff(T_ref)] ;
            
            % interleave the position and joint speed trajectories
            Z_ref = zeros(P.arm_n_states,N_Z) ;
            row_idx = 1 ;
            for idx = 1:2:P.arm_n_states
                Z_ref(idx,:) = Z(row_idx,:) ;
                Z_ref(idx+1,:) = Zd(row_idx,:) ;
                row_idx = row_idx + 1 ;
            end
            
            % generate dummy feedforward output
            U_ref = zeros(P.arm_n_inputs,N_Z) ;
        end
    end
end