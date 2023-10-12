classdef robot_arm_3D_fetch_rnea < robot_arm_rnea_agent
    %% properties
    properties
        % whether or not to plot the CAD version of the arm
        use_CAD_flag = false ;
        link_CAD_data = [] ;
        
        % joint axes
        joint_axes_CAD = [0 0 1 0 1 0 1 1 ;
            0 1 0 1 0 1 0 0 ;
            1 0 0 0 0 0 0 0 ] ;
        
        % joint locations
        joint_locations_CAD = [-0.03 0.12 0.21 0.13 0.20 0.12 0.14 0.16 ;
            0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 ;
            0.75 0.06 0.00 0.00 0.00 0.00 0.00 0.00 ;
            0 0 0 0 0 0 0 0 ;
            0 0 0 0 0 0 0 0 ;
            0 0 0 0 0 0 0 0 ] ;
        
        t_total = 1;
        traj_check_time_discretization = 0.01;
        reference_state = [];
        
        use_mex_controller = true;
    end
    
    %% methods
    methods
        function A = robot_arm_3D_fetch_rnea(varargin)
            use_robust_input = true;
            use_random_interval_params = true;
            
            dimension = 3 ;
                      
            n_links_and_joints = 6 ;
            
            n_states = 2*n_links_and_joints ;
            
            link_shapes = {'cylinder', 'cylinder','cylinder','cylinder','cylinder', 'cylinder'} ;
            
%             link_sizes = [0.05, 0.300, 0.05, 0.200 ;  % size in x
%                           0.05, 0.025, 0.05, 0.025 ;  % size in y
%                           0.05, 0.025, 0.05, 0.025] ; % size in z

            cylinder_diameter = 0.150;
            buffer_dist = cylinder_diameter;
%             buffer_dist = 0.1460*sqrt(2);

            % for rectangular links:
%             link_sizes = [0.1206, 0.4635, 0.001, 0.4254, 0.001, 0.3810 ;  % size in x
%                           0.1460, 0.1460, 0.001, 0.150, 0.001, 0.1460 ;  % size in y
%                           0.0825, 0.1460, 0.001, 0.150, 0.001, 0.1460] ; % size in z
            
            % for cylindrical links:
            % link_sizes = [0.1206, 0.4635, 0.001, 0.4254, 0.001, 0.3810 ;  % size in x
            %               0.1460, 0.1460, 0.001, 0.150, 0.001, 0.1460 ;  % size in y
            %               0.1460, 0.1460, 0.001, 0.150, 0.001, 0.1460] ; % size in z

            % for cylindrical links, updated 20211130 by patrick:
            link_sizes = [0.117, 0.219, 0.133, 0.197, 0.1245, 0.3302 ;  % size in x
                          0.01, cylinder_diameter, cylinder_diameter, cylinder_diameter, cylinder_diameter, cylinder_diameter ;  % size in y
                          0.01, cylinder_diameter, cylinder_diameter, cylinder_diameter, cylinder_diameter, cylinder_diameter] ; % size in z

            joint_state_indices = 1:2:n_states ;
            
            joint_speed_indices = 2:2:n_states ;
            
            joint_types = {'revolute','revolute','revolute','revolute', 'revolute', 'revolute'} ;
            
            joint_axes = [0 0 1 0 1 0;
                          0 1 0 1 0 1;
                          1 0 0 0 0 0] ; % axes are in preceding link's frame
                      
            % joint_locations = [-0.03265, +0.1206/2, +0.3556/2, +0.00, +0.3302/2, +0.00 ; % predecessor x
            %                    +0.00, +0.00, +0.00, +0.00, +0.00, +0.00 ; % predecessor y
            %                    +0.72601, +0.0825/2, +0.00, +0.00, +0.00, +0.00 ; % predecessor z
            %                    -0.1206/2, -0.3556/2, +0.00, -0.3302/2, +0.00, -0.1460 ; % successor x
            %                    +0.00, +0.00, +0.00, +0.00, +0.00, +0.00 ; % successor y
            %                    -0.0825/2, +0.00, +0.00, +0.00, +0.00, +0.00 ];% successor z

            % updated 20211130 by patrick:
            joint_locations = [-0.03260, +0.117/2, +0.219/2, +0.133/2, +0.1970/2, +0.1245/2 ; % predecessor x
                               +0.00000, +0.00000, +0.00000, +0.00000, +0.000000, +0.000000 ; % predecessor y
                               +0.72600, +0.060/2, +0.00000, +0.00000, +0.000000, +0.000000 ; % predecessor z
                               -0.117/2, -0.219/2, -0.133/2, -0.197/2, -0.1245/2, -0.3302/2 ; % successor x
                               +0.00000, +0.00000, +0.00000, +0.00000, +0.000000, +0.000000 ; % successor y
                               -0.060/2, +0.00000, +0.00000, +0.00000, +0.000000, +0.000000 ];% successor z
                           
            joint_state_limits = [-1.6056, -1.221, -Inf, -2.251, -Inf, -2.16;
                                  1.6056, 1.518, Inf, 2.251, Inf, 2.16];
                        
            joint_speed_limits = (Inf).*repmat([-1;1],1,n_links_and_joints) ;
            
%             joint_input_limits = 3.*joint_speed_limits;
            joint_input_limits = [-33.82, -131.76, -76.94, -66.18, -29.35, -25.70;
                                   33.82,  131.76,  76.94,  66.18,  29.35,  25.70];
            
            kinematic_chain = [0 1 2 3 4 5;
                               1 2 3 4 5 6] ;
                           
            gravity_direction = [0;0;-1] ;
            
            set_view_when_animating = true ;
            
            animation_view = 3 ;
            
            % urdf
            robot = importrobot('fetch_arm_reduced.urdf');
            robot.DataFormat = 'col';
            robot.Gravity = [0 0 -9.81];
            
            % primary controller gains
            Kp = 50*eye(6);
            Kd = 10*eye(6);
            Kr = 10*eye(6); % default is 10
%             Kr = 1000*eye(6); % default is 10


            % robust input gains
            % r (which actually isn't r) should be 0.00191 because of \|r\|
            % and K_r
            k_P = 3.2423e+03; % default is 3.2423e+03 lam_min = 0.0017, lam_max = 6.8877, r = 0.0191
            k_I = 100;
            ph_P = 1;
            ph_I = 200;
                        
            Z = [];
            
            % check these with jon!
            ultimate_bound = 0.0191;
            ultimate_bound_position = ultimate_bound/10;
            ultimate_bound_velocity = ultimate_bound*2;
                        
            % get robot params
            [true_params, nominal_params, interval_params] = get_interval_params('fetch_arm_reduced.urdf', use_random_interval_params);            
            
            
            A@robot_arm_rnea_agent('use_robust_input', use_robust_input,...
                'dimension',dimension,'n_links_and_joints',n_links_and_joints,...
                'n_links_and_joints',n_links_and_joints,'n_inputs',n_links_and_joints,...
                'n_states',n_states,...
                'link_shapes',link_shapes,'link_sizes',link_sizes,...
                'joint_state_indices',joint_state_indices,...
                'joint_speed_indices',joint_speed_indices,...
                'joint_types',joint_types,'joint_axes',joint_axes,...
                'joint_locations',joint_locations,...
                'joint_state_limits',joint_state_limits,...
                'joint_speed_limits',joint_speed_limits,...
                'joint_input_limits',joint_input_limits,...
                'kinematic_chain',kinematic_chain,...
                'gravity_direction',gravity_direction,...
                'animation_set_view_flag',set_view_when_animating,...
                'animation_view',animation_view,...
                'buffer_dist', buffer_dist,...
                'animation_time_discretization', 0.01,...
                'true_params', true_params,...
                'nominal_params', nominal_params,...
                'interval_params', interval_params,...
                'Kp', Kp,...
                'Kd', Kd,...
                'Kr', Kr,...
                'k_P', k_P,...
                'k_P', k_P,...
                'k_P', k_P,...
                'k_I', k_I,...
                'ph_P', ph_P,...
                'ph_I', ph_I,...
                'ultimate_bound', ultimate_bound,...
                'ultimate_bound_position', ultimate_bound_position,...
                'ultimate_bound_velocity', ultimate_bound_velocity,...
                'use_robust_input', use_robust_input,...
                'robot', robot,...
                varargin{:}) ;
            
            if A.use_CAD_flag
                A.load_CAD_arm_patch_data()
                A.link_plot_edge_color = [0 0 1] ;
                A.link_plot_edge_opacity = 0 ;
                A.link_plot_face_color = [0.8 0.8 1] ;
                A.link_plot_face_opacity = 1 ;
            end
        end
        
        %% integrator
        function [t_out,z_out] = integrator(A,arm_dyn,t_span,z0)
            dt = 0.01;
            t_span = t_span(1):dt:t_span(2);
            [t_out, z_out] = ode15s(arm_dyn, t_span, [z0; 0]);
            t_out = t_out';
            z_out = z_out(:,1:end-1)';
            
        end

        %% dynamics
        function[q_des, qd_des, qdd_des] =  destraj(A, t, q_0, qd_0, k)
            
            t_plan = A.t_total / 2;
            t_stop = A.t_total - t_plan;
            
            if ~isnan(k)
                if t <= t_plan
                    % compute first half of trajectory
                    q_des = q_0 + qd_0*t + (1/2)*k*t^2;
                    qd_des = qd_0 + k*t;
                    qdd_des = k;
                else
                    % compute trajectory at peak
                    q_peak = q_0 + qd_0*t_plan + (1/2)*k*t_plan^2;
                    qd_peak = qd_0 + k*t_plan;

                    % compute second half of trajectory
                    q_des = q_peak + qd_peak*(t-t_plan) + (1/2)*((0 - qd_peak)/t_stop)*(t-t_plan)^2;
                    qd_des = qd_peak + ((0 - qd_peak)/t_stop)*(t-t_plan);
                    qdd_des = (0 - qd_peak)/t_stop;
                end
            else
                % bring the trajectory to a stop in t_stop seconds
                % trajectory peaks at q_0
                q_peak = q_0;
                qd_peak = qd_0;
                
                if t <= t_stop % we're braking!
                    q_des = q_peak + qd_peak*t + (1/2)*((0 - qd_peak)/t_stop)*t^2;
                    qd_des = qd_peak + ((0 - qd_peak)/t_stop)*t;
                    qdd_des = (0 - qd_peak)/t_stop;
                else % we should already be stopped, maintain that.
                    q_des = q_peak + qd_peak*(t_stop) + (1/2)*((0 - qd_peak)/t_stop)*(t_stop)^2;
                    qd_des = qd_peak + ((0 - qd_peak)/t_stop)*(t_stop); % sanity check, should be 0
                    qdd_des = 0;
                end
            end
        end

        function zd = dynamics(A, t, z, q_0, qd_0, k)    
            % get actual robot state
            z = z(:);
            q = z(A.joint_state_indices);
            qd = z(A.joint_speed_indices);
            int_f = z(end);

            % compute desired trajectory
            [q_des, qd_des, qdd_des] = A.destraj(t, q_0, qd_0, k);

            % error terms
            err = q_des - q;
            d_err = qd_des - qd;

            % modified reference trajectory
            qd_ref = qd_des + A.Kr * err;
            qdd_ref = qdd_des + A.Kr * d_err;

            r = d_err + A.Kr*err;

            % true dynamics 
            [M, C, g] = A.calculate_dynamics(q, qd, A.true_params);
            
            if A.use_mex_controller
%                 fetchMEXController(q, qd, q_des, qd_des, qdd_des, deltaT, eAcc);
                % 7DOF is hard-coded... just set to 0.
                q_mex = [q; 0]; qd_mex = [qd; 0]; q_des_mex = [q_des; 0]; qd_des_mex = [qd_des; 0]; qdd_des_mex = [qdd_des; 0];
                u = fetchMEXController(q_mex, qd_mex, q_des_mex, qd_des_mex, qdd_des_mex, 0, 0);
                u = u(1:end-1);
            else
                % nominal controller
                tau = rnea(q, qd, qd_ref, qdd_ref, true, A.nominal_params);

                % robust input
                if A.use_robust_input
                    v = A.calculate_robust_input(q, qd, qd_ref, qdd_ref, r, tau, int_f);
                else
                    v = zeros(size(tau));
                end
                u_mat = tau + v;
            end

            % update acceleration 
            qdd = M\(u-C*qd-g);

            % preallocate dynamics
            zd = zeros(2*A.n_links_and_joints+1,1) ;

            % compute dynamics
            zd(A.joint_state_indices) = qd(:) ;
            zd(A.joint_speed_indices) = qdd(:) ;

            % update integral action term
            f = norm(err);
            if f < 1e-3
                f = 0;
            end

            % append error
            zd(end) = f; 
        end

        function [M, C, g] = calculate_dynamics(~, q, qd, params)
            M = rnea_mass(q, params);
            C = rnea_coriolis(q, qd, params);
            g = rnea_gravity(q, params);
        end

        function v = calculate_robust_input(A, q, qd, qd_ref, qdd_ref, r, tau, int_f)
            % calculate robust input
            k = A.k_P + A.k_I * int_f;
            ph = A.ph_P + A.ph_I * int_f;

            % calculate max disturbance
            disturbance = rnea(q, qd, qd_ref, qdd_ref, true, A.interval_params) - tau; 

            % bounds on max disturbance
            Phi_min = abs(infimum(disturbance));
            Phi_max = abs(supremum(disturbance));
            norm_rho = norm(max(Phi_min, Phi_max));

            % robust input
            v = (k*norm_rho + ph)*r;
        end 

        function [torques] = compute_torques(A, tout, zout, q_0, qd_0, k)

            torques = zeros(A.n_links_and_joints, length(tout));
            
            for j = 1:length(tout)
                t = tout(j);
                z = zout(:,j);

                % get robot state
                z = z(:);
                q = z(A.joint_state_indices);
                qd = z(A.joint_speed_indices);
                int_f = z(end);

                % compute desired trajectory
                [q_des, qd_des, qdd_des] = A.destraj(t, q_0, qd_0, k);

                % error terms
                err = q_des - q;
                d_err = qd_des - qd;

                % modified reference trajectory
                qd_ref = qd_des + A.Kr * err;
                qdd_ref = qdd_des + A.Kr * d_err;

                r = d_err + A.Kr*err;

                % nominal controller
                tau = rnea(q, qd, qd_ref, qdd_ref, true, A.nominal_params);

                % robust input
                if A.use_robust_input
                    v = A.calculate_robust_input(q, qd, qd_ref, qdd_ref, r, tau, int_f);
                else
                    v = zeros(size(tau));
                end
                u = tau + v;

                torques(:,j) = u;
            end
        end        
        %% move
        function move(A,t_move,T_ref,U_ref,Z_ref,k)
            % method: move(t_move,T_ref,U_ref,Z_ref)
            %
            % Moves the agent for the duration t_move using the nominal
            % inputs U_ref and nominal trajectory Z_ref that are indexed by
            % the nominal time T_ref.
            %
            % This method assumes that the input is zero-order hold, and
            % the input corresponding to the last time index is zero; this
            % is why the last old input is discarded when the input is
            % updated. Similarly, this method assumes that the nominal time
            % starts at 0.
            
            A.vdisp('Moving!',5)
            if A.time(end) == 0
                A.reference_state = A.state;
            end
            
            % set up default reference trajectory
            if nargin < 5
                Z_ref = [] ;
            end
            
            % get the time, input, and reference trajectory to use for
            % moving the agent
            [T_used,U_used,Z_used] = A.move_setup(t_move,T_ref,U_ref,Z_ref) ;
            
            % get the current state
            zcur = A.state(:,end) ;
            q_0 = zcur(A.joint_state_indices);
            qd_0 = zcur(A.joint_speed_indices);
            
            
            switch A.move_mode
                case 'integrator'
                    % call the ode solver to simulate agent
                    [tout,zout] = A.integrator(@(t,z) A.dynamics(t,z,q_0,qd_0,k),...
                                               [0 t_move], zcur) ;

                    [uout] = A.compute_torques(tout, zout, q_0, qd_0, k);
                case 'direct'
                    % don't call the integrator, just assume the agent
                    % perfectly executes the reference trajectory
                    
                    % get the reference trajectory up to time t_move
                    tout = 0:A.integrator_time_discretization:t_move ;
                    [uout,zout] = match_trajectories(tout,T_ref,U_ref,T_ref,Z_ref) ;
                otherwise
                    error('Please set A.move_mode to ''integrator'' or ''direct''!')
            end
            
            A.commit_move_data(tout,zout,T_used,uout,Z_used) ;

        end
        
        function commit_move_data(A,T_state,Z_state,T_used,U_used,Z_used)
            % method: commit_move_data(T_state,Z_state,T_input,U_input,Z_input)
            %
            % After moving the agent, commit the new state and input
            % trajectories, and associated time vectors, to the agent's
            % state, time, input, and input_time properties.
            
            % update the state, time, input, and input time
            A.state = [A.state, Z_state(:,2:end)] ;
            A.time = [A.time, A.time(end) + T_state(2:end)] ;
            A.input_time = [A.input_time, A.input_time(end) + T_used(2:end)] ;
            A.input = [A.input, U_used(:,1:end-1)] ;
            A.reference_state = [A.reference_state, Z_used(:, 2:end)];
        end
        
        %% helper functions   
        function out = input_check(A, t_start)
            % create time vector for checking
            t_agent = A.time(end);
            t_check = t_start:A.traj_check_time_discretization:t_agent;

            if isempty(t_check) || t_check(end) ~= t_agent
                t_check = [t_check, t_agent] ;
            end

            % get agent input trajectory interpolated to time
            u_agent = match_trajectories(t_check,A.time,A.input) ;

            % check torque bounds
            A.vdisp('Running input check!',3);
            out = false;
            for t_idx = 1:length(t_check)
                u = u_agent(:, t_idx);
                for i = 1:length(u)
                    if u(i) > A.joint_input_limits(2, i) || u(i) < A.joint_input_limits(1, i)
                        fprintf('Time %.2f, joint %d torque exceeded: %.2f vs +-%.2f \n', t_check(t_idx), i, u(i), A.joint_input_limits(2, i));
                        out = true;
                    end
                end
            end
            
            if ~out
                A.vdisp('No inputs exceeded', 3);
            end
        end
        
        function out = ultimate_bound_check(A, t_start)
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
                    if abs(q(i) - q_ref(i)) > A.ultimate_bound_position
                        fprintf('Time %.2f, joint %d position bound exceeded: %.5f vs +-%.5f \n', t_check(t_idx), i, abs(q(i) - q_ref(i)), A.ultimate_bound_position);
                        out = true;
                    end
                    if abs(qd(i) - qd_ref(i)) > A.ultimate_bound_velocity
                        fprintf('Time %.2f, joint %d velocity bound exceeded: %.5f vs +-%.5f \n', t_check(t_idx), i, abs(qd(i) - qd_ref(i)), A.ultimate_bound_velocity);
                        out = true;
                    end
                end
            end
            
            if ~out
                A.vdisp('No ultimate bound exceeded', 3);
            end
            
        end
        
        function out = joint_limit_check(A, t_start)
            % create time vector for checking
            t_agent = A.time(end);
            t_check = t_start:A.traj_check_time_discretization:t_agent;

            if isempty(t_check) || t_check(end) ~= t_agent
                t_check = [t_check, t_agent] ;
            end

            % get agent state trajectories interpolated to time
            z_agent = match_trajectories(t_check,A.time,A.state) ;

            % check bound satisfaction
            A.vdisp('Running joint limits check!',3);
            out = false;
            for t_idx = 1:length(t_check)
                q = z_agent(A.joint_state_indices, t_idx);
                qd = z_agent(A.joint_speed_indices, t_idx);
                for i = 1:length(q)
                    if q(i) < A.joint_state_limits(1, i)
                        fprintf('Time %.2f, joint %d position limit exceeded: %.5f vs %.5f \n', t_check(t_idx), i, q(i), A.joint_state_limits(1, i));
                        out = true;
                    end
                    if q(i) > A.joint_state_limits(2, i)
                        fprintf('Time %.2f, joint %d position limit exceeded: %.5f vs %.5f \n', t_check(t_idx), i, q(i), A.joint_state_limits(2, i));
                        out = true;
                    end
                    if qd(i) < A.joint_speed_limits(1, i)
                        fprintf('Time %.2f, joint %d velocity limit exceeded: %.5f vs %.5f \n', t_check(t_idx), i, qd(i), A.joint_speed_limits(1, i));
                        out = true;
                    end
                    if qd(i) > A.joint_speed_limits(2, i)
                        fprintf('Time %.2f, joint %d velocity limit exceeded: %.5f vs %.5f \n', t_check(t_idx), i, qd(i), A.joint_speed_limits(2, i));
                        out = true;
                    end
                end
            end
            
            if ~out
                A.vdisp('No joint limits exceeded', 3);
            end
            
        end
        

        %% plotting setup
        function [faces,vertices] = create_baselink_plot_patch_data(A)
            % create baselink cone for plotting
            [faces,vertices] = make_cuboid_for_patch(0.025, 0.025, 0.025, [0;0;0]) ;
        end
        
        function load_CAD_arm_patch_data(A)
            A.vdisp('Loading CAD files for arm plotting!',1)
            
            shoulder_pan = stlread('shoulder_pan_link_collision.STL') ;
            shoulder_lift = stlread('shoulder_lift_link_collision.STL') ;
            upper_arm = stlread('upperarm_roll_link_collision.STL') ;
            elbow = stlread('elbow_flex_link_collision.STL') ;
            forearm = stlread('forearm_roll_link_collision.STL') ;
            wrist_flex = stlread('wrist_flex_link_collision.STL') ;
            wrist_roll = stlread('wrist_roll_link_collision.STL') ;
            gripper = stlread('gripper_link.STL') ;
            
            temp_link_CAD_data = {shoulder_pan, shoulder_lift,...
                upper_arm, elbow, forearm, wrist_flex, wrist_roll, gripper} ;
            
            % check to make sure the CAD data are patches, not
            % triangulations
            triangulated_flag = false ;
            for idx = 1:8
                current_data = temp_link_CAD_data{idx} ;
                if isa(current_data,'triangulation')
                    triangulated_flag = true ;
                    new_data.faces = current_data.ConnectivityList ;
                    new_data.vertices = current_data.Points ;
                    temp_link_CAD_data{idx} = new_data ;
                end
            end
            
            if triangulated_flag
                A.vdisp('STL read returned a triangulated data format, but we fixed it :)',7)
            end
            
            A.link_CAD_data = temp_link_CAD_data ;
        end
       
        %% plotting
        function plot_links(A,time_or_config)
            if A.use_CAD_flag
                % get the rotations and translations at the current time
                if length(time_or_config) == 1
                    q = match_trajectories(time_or_config,A.time,A.state(A.joint_state_indices,:)) ;
                else
                    q = time_or_config ;
                end
                
                q = [q ; zeros(2,1)] ;
                
                % get transformations for the plot links
                [R,T] = get_link_rotations_and_translations_from_arm_data(q,...
                    A.joint_axes_CAD,A.joint_locations_CAD) ;
                
                % set the number of plot links
                n = 8 ;
                
                % generate plot data for each link
                link_verts = cell(1,n) ;
                for idx = 1:n
                    link_verts{idx} = (R{idx}*A.link_CAD_data{idx}.vertices' + ...
                        T{idx})' ;
                end
                
                if check_if_plot_is_available(A,'links')
                    for idx = 1:n
                        A.plot_data.links(idx).Faces = A.link_CAD_data{idx}.faces ;
                        A.plot_data.links(idx).Vertices = link_verts{idx} ;
                    end
                else
                    link_array = [] ;
                    for idx = 1:n
                        link_data = patch('Faces',A.link_CAD_data{idx}.faces,...
                            'Vertices',link_verts{idx},...
                            'FaceColor',A.link_plot_face_color,...
                            'FaceAlpha',A.link_plot_face_opacity,...
                            'EdgeColor',A.link_plot_edge_color,...
                            'LineWidth',A.link_plot_edge_width,...
                            'EdgeAlpha',A.link_plot_edge_opacity) ;
                        link_array = [link_array, link_data] ;
                    end
                    A.plot_data.links = link_array ;
                    
                    % turn camlight on
                    camlight
                end
            else
                plot_links@robot_arm_rnea_agent(A,time_or_config) ;
            end
        end
        
        %% custom axis limits
        function lims = get_axis_lims(A)
            % figure out the maximum length of the arm
            L = sum(A.link_sizes(1,:)) ;
            
            % create axis limits
            switch A.dimension
                case 2
                    lims = [-L,L,0,L] ;
                case 3
                    %                     switch A.floor_normal_axis
                    %                         case 1
                    %                             lims = [-L, L, -L, L, -L, L] ;
                    %                         case 2
                    %                             lims = [-L, L, 0, L, -L, L] ;
                    %                         case 3
                    %                             lims = [-L,L,-L,L,0,L] ;
                    %                     end
                    lims = [-0.5*L, L, -L, L, -L, 0.75*L] ;
                    
                    % in case base of robot is not a [0;0;0]:
                    lims = lims + [A.joint_locations(1, 1)*ones(1, 2), A.joint_locations(2, 1)*ones(1, 2), A.joint_locations(3, 1)*ones(1, 2)];
                    
                    % make z = 0 the ground
                    lims(5) = 0;
            end
        end
    end
end