classdef uarmtd_agent < robot_arm_agent
    % update 20220126 to use new robust input
    %% properties
    properties
        robot; % rigidBodyTree from robotics toolbox
        params; % stores kinematic and inertial params
        M_min_eigenvalue; % lower bound on eigenvalues of mass matrix
        transmision_inertia;
        q_index; % use this to get index of of joint angle for non-fixed joints
        link_poly_zonotopes; % store link bounding boxes as PZs

        % for plotting CAD:
        use_CAD_flag = false ; % whether or not to plot the CAD version of the arm
        link_CAD_data = [] ;
        joint_locations_CAD;

        % timing
        t_total = 1;
        traj_check_time_discretization = 0.01;
        
        % additional trajectories to log:
        reference_state = [];
        reference_acceleration = [];
        nominal_input = [];
        robust_input = [];
        disturbance = [];
        lyapunov = [];
        r = []; % modified error
        int_lyapunov = [];
        int_disturbance = [];

        % for dealing with measurement noise:
        add_measurement_noise_ = false;
        measurement_noise_pos_;
        measurement_noise_vel_;
        measurement_noise_time_;
        measurement_noise_size_;
    end
    
    %% methods
    methods
        function A = uarmtd_agent(robot, params, varargin)

            robot.DataFormat = 'col';
            robot.Gravity = [0 0 -9.81];
            gravity_direction = [0;0;-1] ;
            
            dimension = 3 ;

            n_links_and_joints = params.nominal.num_joints;
            n_states = 2*params.nominal.num_q;
            q_index = params.nominal.q_index;

            joint_state_indices = 1:2:n_states ;
            joint_speed_indices = 2:2:n_states ;
            joint_types = params.nominal.joint_types';
            joint_axes = params.nominal.joint_axes;

            link_shapes = repmat({'cuboid'}, 1, n_links_and_joints);
            [link_poly_zonotopes, link_sizes, temp_link_CAD_data] = create_pz_bounding_boxes(robot);
                      
            % these joint locations are taken in link-centered body-fixed frames
            % imagine a frame at the center of link 1, and at the center of link 2
            % we want to write the position of joint 2 in link 1's frame, and in link 2's frame.
            for i = 1:n_links_and_joints
                if i == 1
                    joint_locations(:, i) = [params.nominal.T0(1:3, end, i);
                                             -link_poly_zonotopes{i}.c];
                    joint_locations_CAD(:, i) = [params.nominal.T0(1:3, end, i);
                                                 zeros(3, 1)];
                else
                    joint_locations(:, i) = [-link_poly_zonotopes{i-1}.c + params.nominal.T0(1:3, end, i);
                                             -link_poly_zonotopes{i}.c];
                    joint_locations_CAD(:, i) = [params.nominal.T0(1:3, end, i);
                                                 zeros(3, 1)];
                end
            end

            % pull joint position limits from rigidBodyTree
            % MATLAB doesn't store velocity/input limits, so pass these in externally for now.
            for i = 1:n_links_and_joints
                if ~strcmp(robot.Bodies{i}.Joint.Type, 'fixed')
                    joint_state_limits(1, q_index(i)) = robot.Bodies{i}.Joint.PositionLimits(1);
                    joint_state_limits(2, q_index(i)) = robot.Bodies{i}.Joint.PositionLimits(2);
                end
            end
       
            % assuming serial kinematic chain!
            kinematic_chain = [0:n_links_and_joints-1; 1:n_links_and_joints];
            
            set_view_when_animating = true ;
            animation_view = 3 ;

            A@robot_arm_agent('dimension',dimension,'n_links_and_joints',n_links_and_joints,...
                'n_links_and_joints',n_links_and_joints,'n_inputs',n_states/2,...
                'n_states',n_states,...
                'link_shapes',link_shapes,'link_sizes',link_sizes,...
                'joint_state_indices',joint_state_indices,...
                'joint_speed_indices',joint_speed_indices,...
                'joint_types',joint_types,'joint_axes',joint_axes,...
                'joint_locations',joint_locations,...
                'joint_state_limits',joint_state_limits,...
                'kinematic_chain',kinematic_chain,...
                'gravity_direction',gravity_direction,...
                'animation_set_view_flag',set_view_when_animating,...
                'animation_view',animation_view,...
                'animation_time_discretization', 0.01,...
                varargin{:}) ;
            
            A.robot = robot;
            A.params = params;
            A.link_poly_zonotopes = link_poly_zonotopes;
            A.q_index = q_index;
            A.joint_locations_CAD = joint_locations_CAD;
                    
            if A.use_CAD_flag
                A.load_CAD_arm_patch_data(temp_link_CAD_data)
                A.link_plot_edge_color = [0 0 1] ;
                A.link_plot_edge_opacity = 0 ;
                A.link_plot_face_color = [0.8 0.8 1] ;
                A.link_plot_face_opacity = 1 ;
            end
        end
        
        %% info
        function agent_info = get_agent_info(A)
            % call superclass
            agent_info = get_agent_info@robot_arm_agent(A) ;

            agent_info.n_states = A.n_states;
            agent_info.n_inputs = A.n_inputs;
            agent_info.joint_input_limits = A.joint_input_limits;

            % add ref state, inputs, and disturbances
            agent_info.reference_state = A.reference_state ;
            agent_info.reference_acceleration = A.reference_acceleration ;
            agent_info.input = A.input;
            agent_info.nominal_input = A.nominal_input ;
            agent_info.robust_input = A.robust_input ;
            agent_info.disturbance = A.disturbance ;

            % add robot params and LLC info
            agent_info.params = A.params;
            agent_info.link_poly_zonotopes = A.link_poly_zonotopes;
            if isprop(A, 'LLC')
                agent_info.LLC_info = A.LLC.get_LLC_info();
            end
        end

        %% forward kinematics
        function [R,T,J] = get_link_rotations_and_translations(A,time_or_config,cad_flag)
            % [R,T] = A.get_link_rotations_and_translations(time)
            % [R,T] = A.get_link_rotations_and_translations(configuration)
            % [R,T,J] = A.get_link_rotations_and_translations(t_or_q)
            %
            % Compute the rotation and translation of all links in the
            % global (baselink) frame at the given time. If no time is
            % given, then it defaults to 0.
            %
            % The optional third output is the joint locations in 2- or 3-D
            % space, which is also output by A.get_joint_locations(t_or_q).
            %
            % Updated by Patrick on 20220425 to deal with fixed joints
            
            if nargin < 2
                time_or_config = 0 ;
            end
            
            if nargin < 3
                cad_flag = false;
            end
            
            % get joint data
            if size(time_or_config,1) == 1
                t = time_or_config ;
                if t > A.time(end)
                    t = A.time(end) ;
                    warning(['Invalid time entered! Using agent''s final ',...
                        'time t = ',num2str(t),' instead.'])
                end
                
                % interpolate the state for the corresponding time
                z = match_trajectories(t,A.time,A.state) ;
                j_vals = z(1:2:end) ; % joint values
            else
                % assume a configuration was put in
                q = time_or_config ;
                
                if length(q) == A.n_states
                    q = q(1:2:end) ;
                elseif length(q) ~= A.n_states/2
                    error('Please provide either a time or a joint configuration.')
                end
                j_vals = q ;
            end
            
            if ~cad_flag
                j_locs = A.joint_locations ; % joint locations
            else
                j_locs = A.joint_locations_CAD ; % joint locations 
            end
            
            % extract dimensions
            n = A.n_links_and_joints ;
            d = A.dimension ;
            
            % allocate cell arrays for the rotations and translations
            R = mat2cell(repmat(eye(d),1,n),d,d*ones(1,n)) ;
            T = mat2cell(repmat(zeros(d,1),1,n),d,ones(1,n)) ;
            
            % allocate array for the joint locations
            J = nan(d,n) ;
            
            % move through the kinematic chain and get the rotations and
            % translation of each link
            for idx = 1:n
                k_idx = A.kinematic_chain(:,idx) ;
                p_idx = k_idx(1) ;
                s_idx = k_idx(2) ;
                
                % get the rotation and translation of the predecessor and
                % successor links; we assume the baselink is always rotated
                % to an angle of 0 and with a translation of 0
                if p_idx == 0
                    R_pred = eye(d) ;
                    T_pred = zeros(d,1) ;
                else
                    R_pred = R{p_idx} ;
                    T_pred = T{p_idx} ;
                end
                
                % get the location of the current joint
                j_loc = j_locs(:,idx) ;
                
                % compute link rotation
                switch A.joint_types{idx}
                    case 'revolute'
                        % get value of current joint
                        j_idx = j_vals(A.q_index(idx)) ;
                        if d == 3
                            % rotation matrix of current link
%                             axis_pred = A.robot.Bodies{idx}.Joint.JointToParentTransform(1:3, 1:3)*R_pred*A.joint_axes(:,idx) ;
%                             R_succ = axis_angle_to_rotation_matrix_3D([axis_pred', j_idx])*A.robot.Bodies{idx}.Joint.JointToParentTransform(1:3, 1:3)*R_pred ;
                            axis_pred =A.joint_axes(:,idx) ;
                            R_succ = R_pred*A.robot.Bodies{idx}.Joint.JointToParentTransform(1:3, 1:3)*axis_angle_to_rotation_matrix_3D([axis_pred', j_idx]) ;
                        else
                            % rotation matrix of current link
                            R_succ = rotation_matrix_2D(j_idx)*R_pred ;
                        end
                        
                        % create translation
                        T_succ = T_pred + R_pred*j_loc(1:d) - R_succ*j_loc(d+1:end) ;
                    case 'prismatic'
                        % R_succ = R_pred ;
                        error('Prismatic joints are not yet supported!')
                    case 'fixed'
                        if d == 3
                            R_succ = A.robot.Bodies{idx}.Joint.JointToParentTransform(1:3, 1:3)*R_pred ;
                        else
                            % rotation matrix of current link assumed same as predecessor
                            R_succ = R_pred ;
                        end
                        % create translation
                        T_succ = T_pred + R_pred*j_loc(1:d) - R_succ*j_loc(d+1:end) ;
                    otherwise
                        error('Invalid joint type!')
                end
                
                % fill in rotation and translation cells
                R{s_idx} = R_succ ;
                T{s_idx} = T_succ ;
                
                % fill in the joint location
                j_loc_local = j_locs((d+1):end,idx) ;
                J(:,idx) = -R_succ*j_loc_local + T_succ ;
            end
        end
        
        %% integrator
        function [t_out,z_out] = integrator(A,arm_dyn,t_span,z0)
            
            if A.add_measurement_noise_
                A.add_measurement_noise(t_span);
            end

            % setup ODE options
            dt = A.traj_check_time_discretization;
            t_span = t_span(1):dt:t_span(2);
            options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
            
            % call ODE solver
            [t_out, z_out] = ode15s(arm_dyn, t_span, z0, options);
%             [t_out, z_out] = ode45(arm_dyn, t_span, z0, options);

            % process
            t_out = t_out';
            z_out = z_out';            
        end

        % add measurement noise
        function add_measurement_noise(A, t_span)
            A.measurement_noise_time_ = linspace(t_span(1), t_span(end), A.measurement_noise_size_);
            if A.add_measurement_noise_            
%                 A.measurement_noise_pos_ = 5e-4 * randn(A.n_states/2, A.measurement_noise_size_); % noise profile should try to match actual joint encoders
%                 A.measurement_noise_vel_ = 5e-4 * randn(A.n_states/2, A.measurement_noise_size_);
                A.measurement_noise_pos_ = 1e-4 * randn(A.n_states/2, A.measurement_noise_size_); % noise profile should try to match actual joint encoders
                A.measurement_noise_vel_ = 1e-4 * randn(A.n_states/2, A.measurement_noise_size_);
            else            
                A.measurement_noise_pos_ = zeros(A.n_states/2, A.measurement_noise_size_);
                A.measurement_noise_vel_ = zeros(A.n_states/2, A.measurement_noise_size_);
            end        
        end
        
        %% reset
        function reset(A,state,joint_speeds)
            
            A.vdisp('Resetting states',3) ;
            
            % reset to zero by default
            A.state = zeros(A.n_states,1) ;
            A.reference_state = zeros(A.n_states,1) ;
            A.reference_acceleration = zeros(A.n_states/2,1) ;
            
            if nargin > 1
                if length(state) == A.n_states/2
                    % fill in joint positions if they are provided
                    A.vdisp('Using provided joint positions',6)
                    A.state(A.joint_state_indices) = state ;
                    A.reference_state(A.joint_state_indices) = state ;
                    
                    if nargin > 2
                        % fill in joint speeds if they are provided
                        A.vdisp('Using provided joint speeds',6)
                        A.state(A.joint_speed_indices) = joint_speeds ;
                        A.reference_state(A.joint_speed_indices) = joint_speeds ;
                    end
                elseif length(state) == A.n_states
                    % fill in full position and speed state if provided
                    A.vdisp('Using provided full state',6)
                    A.state = state ;
                    A.reference_state = state ;
                else
                    error('Input has incorrect number of states!')
                end
            end
            
            A.vdisp('Resetting time and inputs',3)
            A.time = 0 ;
            A.input = zeros(A.n_inputs,1) ;
            A.input_time = 0 ;
            
            % reset LLC
            if isa(A.LLC,'arm_PID_LLC')
                A.vdisp('Resetting low-level controller integrator error.',3)
                A.LLC.position_error_state = zeros(length(A.joint_state_indices),1) ;
            end
        end       

        %% dynamics
        function zd = dynamics(A, t, z, planner_info)    
            % get actual robot state
            z = z(:);
            q = z(A.joint_state_indices);
            qd = z(A.joint_speed_indices);

            % true dynamics 
            [M, C, g] = A.calculate_dynamics(q, qd, A.params.true); % this is so slow!!!
%             [H, F] = HandC(A.model, q, qd);

            for i = 1:length(q)
                M(i,i) = M(i,i) + A.transmision_inertia(i);
            end
            
            % add measurement noise if desired
            if A.add_measurement_noise_
                [noise_pos, noise_vel] = match_trajectories(t, A.measurement_noise_time_, A.measurement_noise_pos_, A.measurement_noise_time_, A.measurement_noise_vel_, 'linear'); % linearly interpolate noise
                q_meas = q + noise_pos;
                qd_meas = qd + noise_vel;
            else
                q_meas = q;
                qd_meas = qd;
            end
            z_meas = zeros(size(z));
            z_meas(A.joint_state_indices, 1) = q_meas;
            z_meas(A.joint_speed_indices, 1) = qd_meas;

            [u, tau, v] = A.LLC.get_control_inputs(A, t, z_meas, planner_info);

            % update acceleration 
            qdd = M\(u-C*qd-g);
%             qdd = H\(u-F);

            % preallocate dynamics
            zd = zeros(A.n_states,1) ;

            % compute dynamics
            zd(A.joint_state_indices) = qd(:) ;
            zd(A.joint_speed_indices) = qdd(:) ;
        end

        function [M, C, g] = calculate_dynamics(~, q, qd, params)
            M = rnea_mass(q, params);
            C = rnea_coriolis(q, qd, params);
            g = rnea_gravity(q, params);
        end
        
        %% move
        function move(A,t_move,T_ref,U_ref,Z_ref,planner_info)
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
            
            % set up default reference trajectory
            if nargin < 5
                Z_ref = [] ;
            end
            
            % get the time, input, and reference trajectory to use for
            % moving the agent
            [T_used,U_used,Z_used] = A.move_setup(t_move,T_ref,U_ref,Z_ref) ;

            % also get reference acceleration:
            qdd_des = [];
            for i = 1:length(T_used)
                [~, ~, qdd_des(:, i)] = planner_info.desired_trajectory{end}(T_used(i));
            end
            
            % get the current state
            zcur = A.state(:,end) ;
            
            switch A.move_mode
                case 'integrator'
                    % call the ode solver to simulate agent
                    [tout,zout] = A.integrator(@(t,z) A.dynamics(t,z,planner_info),...
                                               [0 t_move], zcur) ;
                                           
                    % initialize trajectories to log
                    uout = zeros(A.n_inputs, size(tout, 2));
                    nominal_out = zeros(size(uout));
                    robust_out = zeros(size(uout));
%                     disturbance_out = zeros(size(uout));
%                     lyap_out = zeros(size(tout));
%                     r_out = zeros(A.n_states/2, size(tout, 2));
%                     int_disturbance_out = interval(zeros(size(uout)), zeros(size(uout)));
%                     int_lyap_out = interval(zeros(size(tout)), zeros(size(tout)));

                    % store approximate inputs at each time:
                    for j = 1:length(tout)
                        t = tout(j);
                        z_meas = zout(:,j);
%                         [uout(:, j), nominal_out(:, j), robust_out(:, j),...
%                          disturbance_out(:, j), lyap_out(:, j), r_out(:, j),...
%                          int_disturbance_out(:, j), int_lyap_out(:, j)] = ...
%                          A.LLC.get_control_inputs(A, t, z_meas, planner_info);
                        [uout(:, j), nominal_out(:, j), robust_out(:, j)] = ...
                         A.LLC.get_control_inputs(A, t, z_meas, planner_info);
                    end
                case 'direct'
                    % don't call the integrator, just assume the agent
                    % perfectly executes the reference trajectory
                    
                    % get the reference trajectory up to time t_move
                    tout = 0:A.integrator_time_discretization:t_move ;
                    zout = match_trajectories(tout,T_ref,Z_ref) ;
                    uout = zeros(A.n_inputs, size(tout, 2));
                    nominal_out = zeros(size(uout));
                    robust_out = zeros(size(uout));
%                     disturbance_out = zeros(size(uout));
%                     lyap_out = zeros(size(tout));
%                     int_disturbance_out = interval(zeros(size(uout)), zeros(size(uout)));
%                     int_lyap_out = interval(zeros(size(tout)), zeros(size(tout)));
%                     r_out = zeros(A.n_states/2, size(tout, 2));
                otherwise
                    error('Please set A.move_mode to ''integrator'' or ''direct''!')
            end
            
%             A.commit_move_data(tout,zout,T_used,uout,Z_used,nominal_out,robust_out,disturbance_out,lyap_out,int_disturbance_out,int_lyap_out,r_out,qdd_des) ;
            A.commit_move_data(tout,zout,T_used,uout,Z_used,nominal_out,robust_out,qdd_des) ;

        end
        
%         function commit_move_data(A,T_state,Z_state,T_used,U_used,Z_used,U_nominal_used,U_robust_used,U_disturbance_used,V_used,int_U_disturbance_used,int_V_used,r_used,qdd_des)
        function commit_move_data(A,T_state,Z_state,T_used,U_used,Z_used,U_nominal_used,U_robust_used,qdd_des)
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
            A.nominal_input = [A.nominal_input, U_nominal_used(:, 1:end-1)];
            A.robust_input = [A.robust_input, U_robust_used(:, 1:end-1)];
%             A.disturbance = [A.disturbance, U_disturbance_used(:, 1:end-1)];
            A.reference_state = [A.reference_state, Z_used(:, 2:end)];
            A.reference_acceleration = [A.reference_acceleration, qdd_des(:, 2:end)];
%             A.lyapunov = [A.lyapunov, V_used(1, 1:end-1)];
%             A.int_disturbance = [A.int_disturbance, int_U_disturbance_used(:, 1:end-1)];
%             A.int_lyapunov = [A.int_lyapunov, int_V_used(1, 1:end-1)];
%             A.r = [A.r, r_used(:, 1:end-1)];
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
        
        function load_CAD_arm_patch_data(A, temp_link_CAD_data)
            A.vdisp('Loading CAD files for arm plotting!',1)
            
            % check to make sure the CAD data are patches, not triangulations
            triangulated_flag = false ;
            for idx = 1:length(temp_link_CAD_data)
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
                                
                % get transformations for the plot links
                cad_flag = true;
                [R,T] = A.get_link_rotations_and_translations(q,cad_flag) ;
                
                % set the number of plot links
                n = A.n_links_and_joints;
                
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
                % show(A.robot, q, 'PreservePlot', false);
                % drawnow;
            else
                plot_links@robot_arm_agent(A,time_or_config) ;
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
                    lims = [-0.5*L, L, -L, L, -L, 0.75*L] ;
                    
                    % in case base of robot is not at [0;0;0]:
                    lims = lims + [A.joint_locations(1, 1)*ones(1, 2), A.joint_locations(2, 1)*ones(1, 2), A.joint_locations(3, 1)*ones(1, 2)];
                    
                    % make z = 0 the ground
                    lims(5) = 0;
            end
        end
    end
end
