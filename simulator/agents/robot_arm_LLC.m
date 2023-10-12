classdef robot_arm_LLC < low_level_controller
    properties
        % agent properties
        arm_dimension
        arm_n_links_and_joints
        arm_joint_state_indices
        arm_joint_speed_indices
        arm_use_robotics_toolbox_model_for_dynamics_flag
    end
    
    methods
        %% constructor
        function LLC = robot_arm_LLC(varargin)
            LLC = parse_args(LLC,varargin{:}) ;
        end
        
        %% setup
        function setup(LLC,agent)
            % call default setup
            setup@low_level_controller(LLC,agent)
            
            % get agent info
            LLC.arm_dimension = agent.dimension ;
            LLC.arm_n_links_and_joints = agent.n_links_and_joints ;
            LLC.arm_joint_state_indices = agent.joint_state_indices ;
            LLC.arm_joint_speed_indices = agent.joint_speed_indices ;
        end
    end
end