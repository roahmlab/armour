classdef robot_arm_RRT_HLP_2D < robot_arm_RRT_HLP
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
         %% constructor
        function HLP = robot_arm_RRT_HLP_2D(varargin)
            
            HLP@robot_arm_RRT_HLP(varargin{:}) ;
            
        end
        
        function plan_path(HLP, start, goal)
            HLP.nodes = [start];
            
        end
        
        %% function generate_waypoints
        
        %% function get_waypoint
        % waypoint = get_waypoint(agent_info,obstacles,lookahead_distance)
        
        %% collision checking
        % obs = vector of obstacles provided by world 
        % arm = Info.collision_check_volume
        
        function out = collision_check_single_obstacle(HLP,obs,arm)
            
            obs = obs.collision_check_patch_data.vertices ;
            
            [x_int,~] = polyxpoly(arm(1,:)',arm(2,:)',...
                obs(:,1),obs(:,2)) ;
            
            if isempty(x_int)
                out = false ;
            else
                out = true ;
            end
        end
        % O = W.obstacles
        function out = collision_check_config(HLP,O,I, q_check)
            % out = collision_check_single_state(W,agent_info,agent_state)
            %
            % Run a collision check for the given state and return true if
            % it is in collision. This gets called by W.collision_check.
            
            %O = W.obstacles ;
            N_O = length(O) ; % in case W.N_obstacles is wrong
            out = false ; % optimism!
            o_idx = 1 ;
            V_arm = I.get_collision_check_volume(q_check) ;
            
            while (o_idx <= N_O) && ~out
                O_idx = O{o_idx} ;
                out = HLP.collision_check_single_obstacle(O_idx,V_arm) ;
                o_idx = o_idx + 1 ;
            end
        end
            
            
      
        
        function check = self_collision_check(HLP, config, arm_volumes)
            % Check self-collisions
            check = 1;
        end
        
        
    end
end

