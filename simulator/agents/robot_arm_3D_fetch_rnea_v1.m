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
    end
    
    %% methods
    methods
        function A = robot_arm_3D_fetch_rnea(varargin)
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

            % robust input gains
            % r (which actually isn't r) should be 0.00191 because of \|r\|
            % and K_r
            k_P = 3.2423e+03; % default is 3.2423e+03 lam_min = 0.0017, lam_max = 6.8877, r = 0.0191
            k_I = 100;
            ph_P = 1;
            ph_I = 200;
            
            Z = [];
            
            use_robust_input = true;
            
            % get robot params
            [true_params, nominal_params, interval_params] = get_interval_params('fetch_arm_reduced.urdf', false);            
            
            
            A@robot_arm_rnea_agent('dimension',dimension,'n_links_and_joints',n_links_and_joints,...
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

        %% dynamics
        function zd = dynamics(A,t,z,T,U,Z)
            A.Z = Z;
            
            % get joint positions
            q = z(A.joint_state_indices);

            % get joint speeds
            qd = z(A.joint_speed_indices) ;
            
            % get integral term
            int_f = z(end);
            
            % enforce speed constraints
            qd = bound_array_elementwise(qd,...
                A.joint_speed_limits(1,:)',A.joint_speed_limits(2,:)') ;
            
            % compute desired trajectory
            destraj = match_trajectories(t, T, Z);          
            q_des = destraj(1:2:12);
            qd_des = destraj(2:2:12);
            qdd_des = zeros(6,1);
            
            % error terms
            err = q_des - q;
            d_err = qd_des - qd;

            % modified reference trajectory
            qd_ref = qd_des + A.Kr * err;
            qdd_ref = qdd_des + A.Kr * d_err;
            
            r = d_err + A.Kr*err;

            % true dynamics 
            [M, C, g] = A.calculate_true_dynamics(q, qd);
  
            % nominal controller
            tau = rnea_passivity(q, qd, qd_ref, qdd_ref, A.nominal_params);

            % robust input
            if A.use_robust_input
                v = calculate_robust_input(A, q, qd, qd_ref, qdd_ref, r, tau, int_f);
                u = tau + v;
            else
                u = tau;
            end
                
            % desired torques
            

            % u = A.LLC.get_control_inputs(A,t,z,T,U,Z) ;
            
            % store control input for debugging
%             A.control_inputs_3 = [A.control_inputs_3 u];
            
            % check torque bounds
%             for i = 1:length(u)
%                if u(i) < A.joint_input_limits(1, i) || u(i) > A.joint_input_limits(2, i)
%                    fprintf('Joint %d torque exceeded: %.3f vs +-%.3f \n', i, u(i), A.joint_input_limits(2, i));
%                end
%             end
          
                     
            % enforce torque constraints
%             u = bound_array_elementwise(u,...
%                 A.joint_input_limits(1,:)',A.joint_input_limits(2,:)') ;
            
            % update acceleration 
            qdd = M\(u-C*qd-g);

            % preallocate dynamics
            zd = zeros(A.n_states+1,1) ;
            
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
        
        %% integrator
        function [t_out,z_out] = integrator(A,arm_dyn,t_span,z0)
            dt = 0.01;
            t_span = t_span(1):dt:t_span(2);
            [t_out, z_out] = ode15s(arm_dyn, t_span, [z0; 0]);
            t_out = t_out';
            z_out = z_out(:,1:end-1)';
            
        end
        
        function u = odefunc(A, t, z)
            T = 0:0.01:1;
            u = zeros(A.n_links_and_joints, length(t));
            
            for i = 1:length(t)
                
                % get joint positions
                q = z(A.joint_state_indices, i);

                % get joint speeds
                qd = z(A.joint_speed_indices, i) ;

                % get integral term
                int_f = z(end, i);

                % enforce speed constraints
                qd = bound_array_elementwise(qd,...
                    A.joint_speed_limits(1,:)',A.joint_speed_limits(2,:)') ;

                % compute desired trajectory
                destraj = match_trajectories(t(:,i), T, A.Z);          
                q_des = destraj(1:2:4);
                qd_des = destraj(2:2:4);
                qdd_des = [0;0];

                % error terms
                err = q_des - q;
                d_err = qd_des - qd;

                % modified reference trajectory
                qd_ref = qd_des + A.Kr * err;
                qdd_ref = qdd_des + A.Kr * d_err;

                r = d_err + A.Kr*err;
                
                % nominal controller
                tau = rnea_passivity(q, qd, qd_ref, qdd_ref, A.nominal_params);

                % robust input
                v = calculate_robust_input(A, q, qd, qd_ref, qdd_ref, r, tau, int_f);

                % store torque inputs
                u(:,i) = tau + v;  
            end
        end
        
        %% move
        function move(A,t_move,T_ref,U_ref,Z_ref)
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
            [T_used,~,Z_used] = A.move_setup(t_move,T_ref,U_ref,Z_ref) ;
            
            % get the current state
            zcur = A.state(:,end) ;
            
            % call the ode solver to simulate agent
            [tout,zout] = A.integrator(@(t,z) A.dynamics(t,z,T_ref,U_ref,Z_ref),...
                                       [0 t_move], zcur) ;
                                   
            uout = A.compute_inputs(tout, zout, T_ref, Z_ref);
            
            A.commit_move_data(tout,zout,T_used,uout,Z_used) ;
        end
        
        %% helper functions
        function [M, C, g] = calculate_true_dynamics(A, q, qd)
            M = rnea_mass(q, A.true_params);
            C = rnea_coriolis(q,qd, A.true_params);
            g = rnea_gravity(q, A.true_params);
        end
        
        function [M, C, g] = calculate_nominal_dynamics(A, q, qd)
            M = rnea_mass(q, A.nominal_params);
            C = rnea_coriolis(q,qd, A.nominal_params);
            g = rnea_gravity(q, A.nominal_params);
        end
        
        function v = calculate_robust_input(A, q, qd, qd_ref, qdd_ref, r, tau, int_f)
            % calculate robust input
            k = A.k_P + A.k_I * int_f;
            ph = A.ph_P + A.ph_I * int_f;
            
            % calculate max disturbance
            disturbance = rnea_passivity(q, qd, qd_ref, qdd_ref, A.interval_params) - tau;
            
            % bounds on max disturbance
            Phi_min = abs(infimum(disturbance));
            Phi_max = abs(supremum(disturbance));
            norm_rho = norm(max(Phi_min, Phi_max));
            
            % robust input
            v = (k*norm_rho + ph)*r;
        end
        
        function input_check = input_check(A)
            disp('hi');
        end
        
        function [torques, exceeds] = compute_inputs(A, tout, zout, T, Z)
            
            torques = zeros(A.num_joints, length(tout));
            
            exceeds = 0;
            
            for j = 1:length(tout)
                t = tout(j);
                z = zout(:,j);
                
                % get robot state
                z = z(:);
                q = z(1:A.num_joints);
                qd = z(A.num_joints+1:2*A.num_joints);
                int_f = z(end);
                
                % compute desired trajectory
%                 [q_des, qd_des, qdd_des] = destraj(t, q_0, qd_0, k, t_total);
                destraj = match_trajectories(t, T, Z);          
                q_des = destraj(1:2:12);
                qd_des = destraj(2:2:12);
                qdd_des = zeros(6,1);

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
                v = calculate_robust_input(q, qd, qd_ref, qdd_ref, r, tau, int_f, A);
                u = tau + v;
                
                % check torque bounds
                for i = 1:A.num_joints
                    if u(i) > A.torque_limits(i)
                        string = strcat('Final Torque exceeded for Joint ', num2str(i));
                        disp(string);
                        exceeds = exceeds + 1;
                    end
                    
                    if u(i) < -A.torque_limits(i)
                        string = strcat('Final Torque exceeded for Joint ', num2str(i));
                        disp(string);
                        exceeds = exceeds + 1;
                    end
                end
                
                torques(:,j) = u;
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