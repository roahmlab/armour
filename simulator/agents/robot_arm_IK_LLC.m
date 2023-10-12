classdef robot_arm_IK_LLC < low_level_controller
    properties

    end
    
    methods
        %% constructor
        function LLC = robot_arm_IK_LLC(varargin)
            LLC@low_level_controller(varargin{:}) ;
        end
        
         %% get control inputs
        function u = get_control_inputs(LLC,~,t,z_cur,varargin)
            % u = LLC.get_control_inputs(agent,t_cur,z_cur,T_ref,U_ref,Z_ref)
            %
            % Given the current time and state, and a reference trajectory
            % as a time and joint locations, compute the control inputs
            % required to track the reference.
            
            % get reference time and feedforward inputs
            T_ref = varargin{1} ;
            U_ref = varargin{2} ;
            
            % get the reference joint location trajectory, which is an
            % arm_dimension * arm_n_joints sized vector
            Z_ref = varargin{3} ;
            
            % convert
        end
    end
end