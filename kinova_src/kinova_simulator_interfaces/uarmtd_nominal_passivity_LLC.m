classdef uarmtd_nominal_passivity_LLC < robot_arm_LLC
    % utilizes nominal params with RNEA
    % adapted from uarmtd_robust_CBF_LLC, just taking out robust part.
    % but keeping ultimate bound checking
    
    properties
    end
    
    methods
        function LLC = uarmtd_nominal_passivity_LLC(varargin)
            LLC = parse_args(LLC,varargin{:}) ;
        end
        
        %% setup
        function setup(LLC,A)
            % call default setup
            setup@robot_arm_LLC(LLC, A);
        end

        %% info
        function info = get_LLC_info(LLC)
            info = struct();
        end

        %% get control inputs
        function [u, tau, v, true_disturbance, true_V, r] = get_control_inputs(LLC, A, t, z_meas, planner_info)
            % u = LLC.get_control_inputs(A, t_cur,z_cur,planner_info)
            %
            % Given the current time and measured state, and a reference trajectory
            % as a function of time within planner_info, compute the control
            % input that attempts to track the given reference.

            % get actual robot state
            z = z_meas(:);
            q = z(A.joint_state_indices);
            qd = z(A.joint_speed_indices);

            % compute desired trajectory
            [q_des, qd_des, qdd_des] = planner_info.desired_trajectory{end}(t);

            % error terms
            err = q_des - q;
            d_err = qd_des - qd;

            % modified reference trajectory
            qd_ref = qd_des + LLC.Kr * err;
            qdd_ref = qdd_des + LLC.Kr * d_err;

            r = d_err + LLC.Kr*err;

            % nominal controller
            tau = rnea(q, qd, qd_ref, qdd_ref, true, A.params.nominal);

            % input is just nominal input, no robust input.
            u = tau;

            % output more for logging if desired
            if nargout > 3
                    disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
                    V_tmp = rnea(q, zeros(A.n_states/2, 1), zeros(A.n_states/2, 1), r, false, A.params.true);
                    V = 0.5*r'*V_tmp;
                    true_disturbance = disturbance;
                    true_V = V;
            end
        end
        
        %% helper functions
        function out = ultimate_bound_check(LLC, A, t_start)
            A.vdisp('No ultimate bound available when not using robust input.', 3)
            out = false;
        end
    end
        
end

