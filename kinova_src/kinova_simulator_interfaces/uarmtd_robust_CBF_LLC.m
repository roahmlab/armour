classdef uarmtd_robust_CBF_LLC < robot_arm_LLC
    % utilizes interval RNEA and CBF arguments to guarantee tracking
    % performance within some ultimate bound
    
    properties
        Kr = 5;
        V_max = 1e-2
        r_norm_threshold = 0;
        alpha_constant = 10;
        alpha = [];
        Kp = [28.1037,28.1037];
        Ki = [4,4];
        maxError = 1e-5;
        use_true_params_for_robust = false;
        use_disturbance_norm = false;
        ultimate_bound;
        ultimate_bound_position;
        ultimate_bound_velocity;
        if_use_mex_controller = false;
        if_use_armour_robust_controller = true;
        saturate_input = false;
    end
    
    methods
        function LLC = uarmtd_robust_CBF_LLC(varargin)
            LLC = parse_args(LLC,varargin{:}) ;
            LLC.alpha = @(h) LLC.alpha_constant*h;
        end
        
        %% setup
        function setup(LLC,A)
            % call default setup
            setup@robot_arm_LLC(LLC, A);

            % using min eigenvalue of robot mass matrix, compute
            % ultimate bound from controller parameters
            if isprop(A, 'M_min_eigenvalue')
                LLC.ultimate_bound = sqrt(2*LLC.V_max/A.M_min_eigenvalue);
                LLC.ultimate_bound_position = (1/LLC.Kr)*LLC.ultimate_bound;
                LLC.ultimate_bound_velocity = 2*LLC.ultimate_bound;
                LLC.vdisp(sprintf('Computed ultimate bound of %.3f', LLC.ultimate_bound), 8);
            else
                warning('No minimum eigenvalue of agent mass matrix specified, can not compute ultimate bound');
            end
        end

        %% info
        function info = get_LLC_info(LLC)
            info.ultimate_bound = LLC.ultimate_bound;
            info.ultimate_bound_position = LLC.ultimate_bound_position;
            info.ultimate_bound_velocity = LLC.ultimate_bound_velocity;

            info.alpha_constant = LLC.alpha_constant;
            info.Kr = LLC.Kr;
        end

        %% get control inputs
        function [u, tau, v, true_disturbance, true_V, r, int_disturbance, int_V] = get_control_inputs(LLC, A, t, z_meas, planner_info)
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

            if ~LLC.if_use_mex_controller
                % nominal controller
                tau = rnea(q, qd, qd_ref, qdd_ref, true, A.params.nominal);
    
                % robust input
                if norm(r) > LLC.r_norm_threshold
                    if LLC.use_true_params_for_robust
                        % calculate true disturbance
                        disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
    
                        % bounds on max disturbance
                        if LLC.use_disturbance_norm
                            norm_rho = norm(disturbance);
                        else
                            rho = abs(r)'*abs(disturbance);
                        end
                        
                        % compute true Lyapunov function
                        V_tmp = rnea(q, zeros(A.n_states/2, 1), zeros(A.n_states/2, 1), r, false, A.params.true);
                        V = 0.5*r'*V_tmp;
                    else
                        % calculate interval disturbance
                        disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.interval) - tau;
                        
                        % bounds on max disturbance
                        if LLC.use_disturbance_norm
                            Phi_min = abs(infimum(disturbance));
                            Phi_max = abs(supremum(disturbance));
                            norm_rho = norm(max(Phi_min, Phi_max));
                        else
                            rho = supremum(abs(r)'*abs(disturbance));
                        end
                        
                        % compute interval Lyapunov function
                        V_tmp = rnea(q, zeros(A.n_states/2, 1), zeros(A.n_states/2, 1), r, false, A.params.interval);
                        V_int = 0.5*r'*V_tmp;
                        V = supremum(V_int);
                    end
                    
                    % assemble input
                    h = -V + LLC.V_max; % value of CBF function
    
                    if LLC.use_disturbance_norm
                        lambda = max(0, -LLC.alpha(h)/norm(r)^2 + norm_rho/norm(r)); % coefficient on r
                    else
                        lambda = max(0, (-LLC.alpha(h) + rho)/norm(r)^2);
                    end
    
                    v = lambda*r; % positive here! because u = tau + v in this code instead of tau - v
                else
                    v = zeros(size(r));
                    V = 0;
                    if ~LLC.use_true_params_for_robust
                        V_int = interval(0, 0);
                        disturbance = interval(0, 0);                    
                    end
                end
    
                % combine nominal and robust inputs
                u = tau + v;
    
                % output more for logging if desired
                if nargout > 3
                    if ~LLC.use_true_params_for_robust
                        true_disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
                        true_V_tmp = rnea(q, zeros(A.n_states/2, 1), zeros(A.n_states/2, 1), r, false, A.params.true);
                        true_V = 0.5*r'*true_V_tmp;
                        int_V = V_int;
                        int_disturbance = disturbance;
                    else
                        if norm(r) == 0
                            disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
                        end
                        true_disturbance = disturbance;
                        true_V = V;
                        int_V = V;
                        int_disturbance = disturbance;
                    end
                end
            else
                model_uncertainty = (A.params.interval.mass_range(2) - A.params.interval.mass_range(1)) / 2;

                if LLC.if_use_armour_robust_controller
                    [u, tau, v] = kinova_controller(LLC.Kr * ones(A.n_inputs,1), LLC.alpha_constant, LLC.V_max, LLC.r_norm_threshold, q, qd, q_des, qd_des, qdd_des, model_uncertainty);
                else
                    [u, tau, v] = kinova_controller_ALTHOFF(LLC.Kr * ones(A.n_inputs,1), LLC.Kp, LLC.Ki, LLC.maxError, q, qd, q_des, qd_des, qdd_des, model_uncertainty);
                end

                if nargout > 3
                    if ~LLC.use_true_params_for_robust
                        true_disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
                        true_V_tmp = rnea(q, zeros(A.n_states/2, 1), zeros(A.n_states/2, 1), r, false, A.params.true);
                        true_V = 0.5*r'*true_V_tmp;
                    else
                        if norm(r) == 0
                            disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.params.true) - tau;
                        end
                        true_disturbance = disturbance;
                        true_V = V;
                    end
                end
            end

            if LLC.saturate_input
                u = min(max(u, A.joint_input_limits(1,:)'), A.joint_input_limits(2,:)');
            end
        end
        
        %% helper functions
        function out = ultimate_bound_check(LLC, A, t_start)
            % create time vector for checking
            t_agent = A.time(end);
            t_check = t_start:A.traj_check_time_discretization:t_agent;

            if isempty(t_check) || t_check(end) ~= t_agent
                t_check = [t_check, t_agent] ;
            end

            % get agent state and reference trajectories interpolated to time
            z_agent = match_trajectories(t_check,A.time,A.state) ;
            z_ref_agent = match_trajectories(t_check,A.time,A.reference_state) ;

            % check bound satisfaction
            A.vdisp('Running ultimate bound check!',3);
            out = false;
            for t_idx = 1:length(t_check)
                q = z_agent(A.joint_state_indices, t_idx);
                qd = z_agent(A.joint_speed_indices, t_idx);
                q_ref = z_ref_agent(A.joint_state_indices, t_idx);
                qd_ref = z_ref_agent(A.joint_speed_indices, t_idx);
                for i = 1:length(q)
                    if abs(q(i) - q_ref(i)) > LLC.ultimate_bound_position
                        fprintf('        LLC: Time %.2f, joint %d position bound exceeded: %.5f vs +-%.5f \n', t_check(t_idx), i, q(i) - q_ref(i), LLC.ultimate_bound_position);
                        out = true;
                    end
                    if abs(qd(i) - qd_ref(i)) > LLC.ultimate_bound_velocity
                        fprintf('        LLC: Time %.2f, joint %d velocity bound exceeded: %.5f vs +-%.5f \n', t_check(t_idx), i, qd(i) - qd_ref(i), LLC.ultimate_bound_velocity);
                        out = true;
                    end
                end
            end

            if ~out
                LLC.vdisp('No ultimate bound exceeded', 3);
            end
        end
    end
        
end

