classdef robot_arm_PID_LLC < robot_arm_LLC
    properties
        % gains
        K_ff
        K_p
        K_d
        K_i
        
        % method type for interpolation - these are the same as the
        % interpolation methods for MATLAB's interp1, available here:
        %   https://www.mathworks.com/help/matlab/ref/interp1.html#btwp6lt-1-method
        %
        % for zero-order hold, use 'previous', and for first-order hold,
        % use 'linear'
        interp_method = 'linear' ;
        
        % position error for integrator term
        position_error_state
    end
    
    methods
        %% constructor
        function LLC = robot_arm_PID_LLC(varargin)
            LLC@robot_arm_LLC(varargin{:}) ;
        end
        
        %% setup
        function setup(LLC,agent)
            % call default setup
            setup@robot_arm_LLC(LLC,agent)
            
            % create default control gains, assuming the agent's states are
            % in the format
            %   (q_1, qdot_1, q_2, qdot_2, ... , q_n, qdot_n),
            % where q_i is the configuration of joint i
            
            if isempty(LLC.K_ff)
                LLC.K_ff = LLC.make_ff_gain_array(1) ;
            end
            
            if isempty(LLC.K_p)
                LLC.K_p = LLC.make_p_gain_array(100) ;
            end
            
            if isempty(LLC.K_d)
                LLC.K_d = LLC.make_d_gain_array(10) ;
            end
            
            if isempty(LLC.K_i)
                LLC.K_i = LLC.make_i_gain_array(0.01) ;
            end
            
            % reset the integrator error
            LLC.position_error_state = zeros(length(LLC.arm_joint_state_indices),1) ;
        end
        
        %% get control inputs
        function u = get_control_inputs(LLC,~,t,z_cur,varargin)
            % u = LLC.get_control_inputs(~,t_cur,z_cur,T_ref,U_ref,Z_ref)
            %
            % Given the current time and state, and a reference trajectory
            % as a time, input, and state arrays, compute the PID control
            % input that attempts to track the given reference.
            
            % get reference time and (feedforward) input
            T_ref = varargin{1} ;
            U_ref = varargin{2} ;
            
            % get reference input and state at the current time
            if length(varargin) > 2 && ~isempty(varargin{3})
                Z_ref = varargin{3} ;
                [u_ref,z_ref] = match_trajectories(t,T_ref,U_ref,T_ref,Z_ref,LLC.interp_method) ;
            else
                u_ref = match_trajectories(t,T_ref,U_ref) ;
                z_ref = z_cur ;
            end
            
            % compute control input
            pd_error = z_cur - z_ref ;
            i_error = LLC.position_error_state ;
            u = LLC.K_ff*u_ref - LLC.K_p*pd_error + ...
                - LLC.K_d*pd_error - LLC.K_i*i_error ;
            
            
            % add gravity compensation
            if LLC.arm_use_robotics_toolbox_model_for_dynamics_flag
                q = z_cur(LLC.arm_joint_state_indices) ;
                u = u + A.robotics_toolbox_model.gravityTorque(q) ;
            end
            
            % update integrator error
            LLC.position_error_state = i_error + pd_error(LLC.arm_joint_state_indices) ;
        end
        
        %% set control gains
        function set_gains(LLC,varargin)
            % LLC.set_gains(K_ff,K_p,K_i,K_d)
            % LLC.set_gains('FF',K_ff,'P',K_p,'D',K_d,'I',K_i)
            %
            % Set the controller's gains.
            
            % if the first argument in is a string, then we'll parse the
            % arguments as an options list, otherwise we assume that the
            % gains are entered in the order K_ff, K_p, K_i, K_d
            if ischar(varargin{1})
                for idx = 1:2:length(varargin)
                    switch lower(varargin{idx})
                        case 'ff'
                            LLC.K_ff = LLC.make_ff_gain_array(varargin{idx+1}) ;
                        case 'p'
                            LLC.K_p = LLC.make_p_gain_array(varargin{idx+1}) ;
                        case 'i'
                            LLC.K_i = LLC.make_i_gain_array(varargin{idx+1}) ;
                        case 'd'
                            LLC.K_d = LLC.make_d_gain_array(varargin{idx+1}) ;
                        otherwise
                            error('Please pick FF, P, I, or D as the keyword arguments.')
                    end
                end
            else
                LLC.K_ff = LLC.make_ff_gain_array(varargin{1}) ;
                
                if length(varargin) > 1
                    LLC.K_p = LLC.make_p_gain_array(varargin{2}) ;
                end
                
                if length(varargin) > 2
                    LLC.K_i = LLC.make_i_gain_array(varargin{3}) ;
                end
                
                if length(varargin) > 3
                    LLC.K_d = LLC.make_d_gain_array(varargin{4}) ;
                end
            end
        end
        
        function K_ff = make_ff_gain_array(LLC,k_ff)
            if length(k_ff) == 1
                K_ff = k_ff.*eye(LLC.n_agent_inputs) ;
            else
                K_ff = diag(k_ff) ;
            end
        end
        
        function K_p = make_p_gain_array(LLC,k_p)
            if length(k_p) == 1
                K_p = k_p*eye(LLC.n_agent_inputs) ;
            else
                K_p = diag(k_p) ;
            end
            K_p = [K_p ; zeros(LLC.n_agent_inputs)] ;
            K_p = reshape(K_p,LLC.n_agent_inputs,[]) ;
        end
        
        function K_i = make_i_gain_array(LLC,k_i)
            if length(k_i) == 1
                K_i = k_i*eye(LLC.n_agent_inputs) ;
            else
                K_i = diag(k_i) ;
            end
        end
        
        function K_d = make_d_gain_array(LLC,k_d)
            if length(k_d) == 1
                K_d = k_d*eye(LLC.n_agent_inputs) ;
            else
                K_d = diag(k_d) ;
            end
            K_d = [zeros(LLC.n_agent_inputs) ; K_d] ;
            K_d = reshape(K_d,LLC.n_agent_inputs,[]) ;
        end
    end
end