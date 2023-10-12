classdef robot_arm_rnea_agent < multi_link_agent
    properties
        %% DH parameters
        dh_table
        
        %% robot params
        true_params
        nominal_params
        interval_params
        pz_interval_params
        pz_nominal_params
  
        %% controller gains
        % primary controller gains
        Kp
        Kd
        Kr

        % robust input gains
        k_P
        k_I
        ph_P
        ph_I
        
        % new robust input params
        c
        alpha
        alpha_c
        
        % ultimate bound
        ultimate_bound
        ultimate_bound_position
        ultimate_bound_velocity
        
        use_robust_input
        
        % error
        error_indices
        
        % control values
%         control_inputs
%         control_inputs_2
%         control_inputs_3
        Z
        reference_state
        
        %% links
        % the number of links and joints of the arm; we use the convention
        % that every joint is succeeded by a link (and every link except
        % the baselink is preceded by a joint), hence the number of joints
        % and links are the same
        n_links_and_joints
        
        % pick link shapes as a cell array containing the strings:
        % 'box' (2D), 'oval' (2D), 'cuboid' (3D), or ellipsoid (3D)
        % (note that the baselink is plotted as a triangle in 2D and a cone
        % in 3D by default)
        link_shapes
        
        % link sizes is a 2-by-n_links_and_joints array, where each row
        % specifies the breadth of the link along the x, y, and z
        % dimensions; note that, for ellipsoid-shaped links, these values
        % specify the diameters of the ellipsoid along the link's local x,
        % y, and z dimensions; for cylindrical links, the x value is the
        % length of the cylinder and the y and z values are the diameter
        % (only the y value is used to generate volumes)
        link_sizes
        
        % the links all have 1 kg mass by default
        link_masses
        
        % for a 3-D arm, each link has a 3-by-3 inertia matrix stored in a
        % 1-by-n_links_and_joints cell array
        link_inertia_matrices
        
        %% joints
        % note that each link can only be associated with one joint
        % preceding it in the current forward kinematics formulation
        
        % joint types are a cell array containing the strings 'revolute'
        % and 'prismatic' of size 1-by-n_links_and_joints; note
        joint_types
        
        % the joint axes are in the predecessor's coordinate frame; for a
        % 2-D arm, revolute joints should use the z-axis, which is the
        % default here
        joint_axes
        
        % specify the location of each joint as an (x,y) coordinate in the
        % local frame of the pred/succ links; this array must have as many
        % columns as there are joints; note that, for any joint that is
        % succeeded by a link with no successors, the "successor" joint
        % location is the end effector location on that link
        joint_locations
        
        % joint limits and speed limits, where each column corresponds to
        % one joint; the top row is the minimum and the bottom is the
        % maximum, in rad or rad/s respectively
        joint_state_limits
        joint_speed_limits
        
        % state index data
        joint_state_indices
        joint_speed_indices
                
        % specify the kinematic chain of predecessor and successor links;
        % each column corresponds to a joint, the first row is the joint
        % predecessor link, and the second row is the successor link; 0
        % indicates the baselink
        kinematic_chain
        
        %% physics
        % whether or not to simulate gravity
        gravity_flag = false ; % not implemented yet
        gravity_direction
        
        % torque limits; each column corresponds to one joint, the first
        % row is the minimum, and the second row is the maximum in N*m
        joint_input_limits
        
        % amount of time to execute a stopping maneuver
        t_stop = 0.5 ;
        
        %% robotics toolbox model
        % the agent file will autogenerate a robotics toolbox robot model
        % that can be used for the forward dynamics
        robotics_toolbox_model
        use_robotics_toolbox_model_for_dynamics_flag = false ;
        robot;
        %% miscellaneous
        % move mode
        move_mode = 'integrator' ; % choose 'direct' or 'integrator'
        integrator_time_discretization = 0.01 ; % s
        
        % collision checking
        collision_check_patch_data
        
        % plotting
        link_plot_data
        link_plot_face_color = [0 0 1] ;
        link_plot_face_opacity = 0.1 ;
        link_plot_edge_color = [0 0 1] ;
        link_plot_edge_opacity = 1.0 ;
        link_plot_edge_width = 1.25 ;
        
        % floor orientation
%         floor_normal_axis = 3 ;
        
        % buffer distance for obstacles
        buffer_dist = 0 ;
    end
    
    methods
        %% constructor
        function A = robot_arm_rnea_agent(varargin)
            % A = robot_arm_rnea_agent('property1',value1,'property2',value2,...)
            %
            % This class implements a generic multilink robot arm. By
            % default, it is a 2-D, 2-link, 2-DOF arm.
            
            % make sure there are input arguments
            if nargin == 0
                error(['Please create a subclass of this class with all ',...
                    'required properties defined!'])
            end
            
            % by default, the arm's states are (position,speed) of each
            % joint, in order from the first joint onwards
            n_states = 2*A.n_links_and_joints ;
            joint_state_indices = 1:2:n_states ;
            joint_speed_indices = 2:2:n_states ;
            
            n_inputs = A.n_links_and_joints ; % inputs are joint torques by default
            
            % default low-level controller
            LLC = robot_arm_PID_LLC() ;
            
            % create agent
            A = parse_args(A,'n_states',n_states,'n_inputs',n_inputs,...
                'joint_state_indices',joint_state_indices,...
                'joint_speed_indices',joint_speed_indices,...
                'LLC',LLC,...
                'animation_set_axes_flag',true,...
                'animation_gif_filename','arm_animation.gif',...
                'position_indices',[],varargin{:}) ;
            
            % ensure properties are consistent
            A.check_and_fix_properties() ;
            
            % set time, states, and inputs to zero
            A.reset() ;
            
            % set up low-level controller
            A.LLC.setup(A)
            
            % set up plotting and collision check data
            A.create_plot_patch_data() ;
            A.create_collision_check_patch_data() ;
            A.plot_data.baselink = [] ;
            A.plot_data.links = [] ;
        end
        
        %% property check
        function check_and_fix_properties(A)
            % A.check_and_fix_properties()
            %
            % Go through each of the arm's default properties and make sure
            % they are consisten with the total number of links and joints,
            % and with the dimension of the arm's workspace.
            
            A.vdisp('Ensuring all arm properties are consistent',3)
            
            % get the actual number of links and joints; the number of
            % links MUST be equal to the number of joints, because every
            % joint is connected to a single predecessor and a single
            % successor link
            N = A.n_links_and_joints ;
            if isempty(N)
                A.vdisp('Setting A.n_links_and_joints based on link sizes',1)
                N = size(A.link_sizes,2) ;
                A.n_links_and_joints = N ;
            end
            
            % check dimension
            d = A.dimension ;
            if (d ~= 2) && (d ~= 3)
                error(['The arm''s workspace dimension (',num2str(d),...
                    ') is incorrect! Pick 2 or 3.'])
            end
            
            %% states and state indices
            A.vdisp('Checking states and state indices',5)
            
            % by default, the arm's states are (position,speed) of each
            % joint
            if isempty(A.n_states)
                A.n_states = 2*N ;
                A.vdisp(['Setting number of states to ',num2str(A.n_states)],6)
            end
            
            if isempty(A.joint_state_indices)
                A.joint_state_indices = 1:2:A.n_states ;
                A.vdisp('Setting default joint state indices',6)
            end
            
            if isempty(A.joint_speed_indices)
                A.joint_speed_indices = 2:2:A.n_states ;
                A.vdisp('Setting default joint speed indices',6)
            end
            
            if isempty(A.n_inputs)
                A.n_inputs = N ;
                A.vdisp(['Setting number of inputs to ',num2str(A.n_states)],6)
            end
            
            %% links
            A.vdisp('Checking links',5)
            
            % check dimension and number of links
            [d_l, N_l] = size(A.link_sizes) ;
            
            if N_l ~= N
                error(['The arm does not have the same number of links as ',...
                    'its n_links_and_joints property!'])
            end
            
            if d_l ~= d
                error(['The arm''s links are not of the same dimension (',...
                    num2str(d_l),' as the arm''s dimension property (',...
                    num2str(d),')!'])
            end
            
            % check link shapes
            N_link_shapes = length(A.link_shapes) ;
            if N_link_shapes > N
                A.vdisp('Removing extra link shapes',6)
                A.link_shapes = A.link_shapes(1:N) ;
            elseif N_link_shapes < N
                A.vdisp('Setting missing link shapes',6)
                switch d
                    case 2
                        A.link_shapes = [A.link_shapes, repmat({'box'},1,N-N_link_shapes)] ;
                    case 3
                        A.link_shapes = [A.link_shapes, repmat({'cuboid'},1,N-N_link_shapes)] ;
                end
            end
            
            % check link masses
            if isempty(A.link_masses) || length(A.link_masses) ~= N
                A.vdisp('Setting link masses',6)
                A.link_masses = ones(1,N) ;
            end
            
            % check link inertia matrices
            % see: en.wikipedia.org/wiki/List_of_moments_of_inertia
            if isempty(A.link_inertia_matrices)
                A.vdisp('Creating link inertia matrices!',6)
                
                J_cell = cell(1,N) ;
                
                for idx = 1:N
                    m = A.link_masses(idx) ;
                    l = A.link_sizes(:,idx) ;                    
                    switch A.link_shapes{idx}
                        case {'box', 'oval'}
                            % treat things as thin boxes in 2-D
                            I_x = 0 ;
                            I_y = 0 ;
                            I_z = (1/12)*m*sum(l.^2) ;                            
                        case 'cuboid'
                            I_x = (1/12)*m*((l(2)+l(3)).^2) ;
                            I_y = (1/12)*m*((l(1)+l(3)).^2) ;
                            I_z = (1/12)*m*((l(1)+l(2)).^2) ;                            
                        case 'ellipsoid'
                            l = l./2 ;
                            I_x = (1/5)*m*((l(2)+l(3)).^2) ;
                            I_y = (1/5)*m*((l(1)+l(3)).^2) ;
                            I_z = (1/5)*m*((l(1)+l(2)).^2) ;
                        case 'cylinder'
                            if l(2) ~= l(3)
                                warning(['The size definition of link ',num2str(idx),...
                                    ' has two mismatched diameters!'])
                            end
                            
                            h = l(1) ;
                            r = l(2)/2 ;
                            
                            I_x = (1/2)*m*r^2 ;
                            I_y = (1/12)*m*(3*r^2 + h) ;
                            I_z = (1/12)*m*(3*r^2 + h) ;
                        otherwise
                            error([A.link_shapes{idx}, 'is an unsupported link shape!'])
                    end
                    J_cell{idx} = diag([I_x, I_y, I_z]) ;                    
                end                
                A.link_inertia_matrices = J_cell ;
            end
            
            %% joints
            A.vdisp('Checking joints',5)
            
            % check dimension and number of joints
            [d_j_2x, N_j] = size(A.joint_locations) ;
            d_j = d_j_2x / 2 ;
            
            if N_j ~= N
                error(['The arm does not have the same number of joints as ',...
                    'its n_links_and_joints property!'])
            end

            if d_j ~= d
                error(['The arm''s joints are not of the same dimension (',...
                    num2str(d_j),') as the arm''s dimension property (',...
                    num2str(d),')!'])
            end
            
            % set default joint types (all revolute)
            if isempty(A.joint_types) || (length(A.joint_types) ~= N)
                A.joint_types = repmat({'revolute'},1,N) ;
            end
            
            % set default joint axes (note that these axes are in the
            % local coordinate frame of the joint's predecessor link)
            N_joint_axes = size(A.joint_axes,2) ;
            if N_joint_axes > N
                A.vdisp('Removing extra joint axes',6)
                A.joint_axes = A.joint_axes(:,1:N) ;
            elseif N_joint_axes < N
                A.vdisp('Setting missing joint axes',6)
                A.joint_axes = [A.joint_axes, repmat([0;0;1],1,N-N_joint_axes)] ;
            end
            
            if isempty(A.kinematic_chain)
                A.vdisp('Setting arm kinematic chain',6)
                A.kinematic_chain = [0:(N-1) ; 1:N] ;
            end
            
            %% physics
            if isempty(A.gravity_direction)
                A.vdisp('Setting default gravity direction',9)
                switch d
                    case 2
                        A.gravity_direction = [0;-1;0] ;
                    case 3
                        A.gravity_direction = [0;0;-1] ;
                end
            end
            
            %% robotics toolbox model
            if isempty(A.robotics_toolbox_model)
                A.vdisp('Making Robotics Toolbox model',5)
                create_robotics_toolbox_model(A) ;
            end
        end
        
        %% random state generation
        function varargout = create_random_state(A)
            % [z,u] = A.create_random_state()
            % [q,qd,u] = A.create_random_state()
            %
            % Create a random state z and random input u, given the arm's
            % state, speed, and input limits. If three output args are
            % specified, then return the config q, joint speeds qd, and
            % input u.
            
            % set any joint limits that are +Inf to pi and -Inf to -pi
            state_lims = A.joint_state_limits ;
            joint_limit_infs = isinf(state_lims) ;
            state_lims(1,joint_limit_infs(1,:)) = -pi ;
            state_lims(2,joint_limit_infs(2,:)) = +pi ;
            
            % make random state
            q = rand_range(state_lims(1,:),state_lims(2,:))' ;
            qd = rand_range(A.joint_speed_limits(1,:),A.joint_speed_limits(2,:))' ;
            
            % make random input
            u = rand_range(A.joint_input_limits(1,:),A.joint_input_limits(2,:))' ;
            
            switch nargout
                case 1
                    varargout = {q} ;
                case 2                    
                    z = nan(A.n_states,1) ;
                    z(A.joint_state_indices) = q ;
                    z(A.joint_speed_indices) = qd ;
                    
                    varargout = {z, u} ;
                case 3
                    varargout = {q, qd, u} ;
            end
        end
        
        %% plot data
        function create_plot_patch_data(A)
            % A.create_plot_patch_data()
            %
            % This method fills in A.link_plot_data with a cell array of
            % faces vectors (each 1-by-NF) and a cell array of vertices
            % arrays (each NV-by-2).
            
            % set up cell array to save patch data
            plot_faces_cell = cell(1,A.n_links_and_joints) ;
            plot_verts_cell = cell(1,A.n_links_and_joints) ;
            
            % get link sizes
            L = A.link_sizes ;
            
            % make baselink plot data
            [baselink_faces,baselink_vertices] = A.create_baselink_plot_patch_data() ;
                    
            % create links
            for l_idx = 1:size(L,2)
                l = L(:,l_idx) ;
                
                % create link based on link type
                switch A.link_shapes{l_idx}
                    case 'box'
                        [link_faces, link_vertices] = make_box(l) ;
                    case 'oval'
                        [link_faces, link_vertices] = make_oval(l) ;
                    case 'cuboid'
                        % create box for link that is slightly shorter than
                        % the actual volume of the box for prettiness
                        % purposes
                        % l(1) = l(1) - 0.045 ;
                        [link_faces,link_vertices] = make_cuboid_for_patch(l) ;
                    case 'ellipsoid'
                        % remember that the link sizes specify the
                        % diameters of the link in each dimension
                        l = l./2 ;
                        [link_faces,link_vertices] = make_ellipsoid_for_patch(l(1),l(2),l(3),zeros(3,1),6) ;
                    case 'cylinder'
                        % l = (length, radius, (not used))
                        [link_faces,link_vertices] = make_cylinder_for_patch(l(2)/2,l(1),10,true,true) ;
                        R = axang2rotm([0 1 0 pi/2]) ;
                        link_vertices = (R*link_vertices')' ;
                    otherwise
                        error('Invalid link type! Pick cuboid, ellipsoid, or cylinder.')
                end

                % fill in cell array
                plot_faces_cell{l_idx} = link_faces ;
                plot_verts_cell{l_idx} = link_vertices ;
            end
            
            % fill in plot link data object
            A.link_plot_data.link_faces = plot_faces_cell ;
            A.link_plot_data.link_vertices = plot_verts_cell ;
            A.link_plot_data.baselink_faces = baselink_faces ;
            A.link_plot_data.baselink_vertices = baselink_vertices ;
        end
        
        function [faces,vertices] = create_baselink_plot_patch_data(A)
            switch A.dimension
                case 2
%                     % create baselink triangle for plotting
%                     vertices = 0.05.*[-1, 1, 0 ;
%                         0, 0, 1 ]' ;
%                     faces = [1 2 3 1] ;
                    % create tiny square for plotting
                    vertices = 0.025.*[-1, 1, 1, -1 ;
                        -1, -1, 1, 1]' ;
                    faces = [1 2 3 4 1] ;
                    
                case 3
                    % create baselink cone for plotting
                    [faces,vertices] = make_cone_for_patch(0.05,0.05) ;
            end
        end
        
        %% collision check data
        function create_collision_check_patch_data(A)
            % A.create_collision_check_patch_data()
            %
            % Create a cell array of faces vectors (each 1-by-NF) and a
            % cell array of vertices arrays (each NV-by-2) to be used for
            % collision checking.
            
            
            switch A.dimension
                case 2
                    cc_faces_cell = A.link_plot_data.link_faces ;
                    cc_verts_cell = A.link_plot_data.link_vertices ;
                case 3
                    % set up cell array to save patch data
                    cc_faces_cell = cell(1,A.n_links_and_joints) ;
                    cc_verts_cell = cell(1,A.n_links_and_joints) ;
                    
                    % create links as rectangles
                    L = A.link_sizes ;
                    
                    for l_idx = 1:size(L,2)
                        l = L(:,l_idx) ;
                        
                        switch A.link_shapes{l_idx}
                            case 'cuboid'
                                [~,link_vertices] = make_cuboid_for_patch(l) ;
                            case 'ellipsoid'
                                l = l./2 ;
                                [~,link_vertices] = make_ellipsoid_for_patch(l(1),l(2),l(3),zeros(3,1),6) ;
                            case 'cylinder'
                                % l = (length, radius, (not used))
                                [~,link_vertices] = make_cylinder_for_patch(l(2)/2,l(1),10,true,true) ;
                                R = axang2rotm([0 1 0 pi/2]) ;
                                link_vertices = (R*link_vertices')' ;
                            otherwise
                                error('Invalid link type! Pick cuboid, ellipsoid, or cylinder!')
                        end
                        link_faces = convhull(link_vertices) ;
                        
                        % fill in cell array
                        cc_faces_cell{l_idx} = link_faces ;
                        cc_verts_cell{l_idx} = link_vertices ;
                    end
            end
            
            % fill in collision check data object
            A.collision_check_patch_data.faces = cc_faces_cell ;
            A.collision_check_patch_data.vertices = cc_verts_cell ;
        end
        
        %% reset
        function reset(A,state,joint_speeds)
            
            A.vdisp('Resetting states',3) ;
            
            % reset to zero by default
            A.state = zeros(A.n_states,1) ;
            A.reference_state = zeros(A.n_states,1) ;
            
            if nargin > 1
                if length(state) == A.n_links_and_joints
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
        
        %% get agent info
        function agent_info = get_agent_info(A)
            % call superclass
            agent_info = get_agent_info@agent(A) ;

            % add position info
            agent_info.position = A.state(A.joint_state_indices,:) ;
            agent_info.position_indices = A.joint_state_indices ;
            
            % properties
            agent_info.n_links_and_joints = A.n_links_and_joints ;
            agent_info.dimension = A.dimension ;
            agent_info.joint_axes = A.joint_axes ;
            agent_info.joint_locations = A.joint_locations ;
            agent_info.joint_state_indices = A.joint_state_indices ;
            agent_info.joint_speed_indices = A.joint_speed_indices ;
            agent_info.joint_state_limits = A.joint_state_limits ;
            agent_info.joint_speed_limits = A.joint_speed_limits ;
            agent_info.joint_input_limits = A.joint_input_limits ;
            agent_info.reach_limits = A.get_axis_lims() ;
            agent_info.buffer_dist = A.buffer_dist ;
            agent_info.dh_table = A.dh_table;

            % collision check data
            agent_info.collision_check_patch_data = A.collision_check_patch_data ;
            
            % useful functions
            agent_info.get_collision_check_volume = @(q) A.get_collision_check_volume(q) ;
            agent_info.get_joint_locations = @(q) A.get_joint_locations(q) ;
            agent_info.get_link_rotations_and_translations = @(t_or_q) A.get_link_rotations_and_translations(t_or_q) ;
            agent_info.inverse_kinematics = @(J,q_0) A.inverse_kinematics(J,q_0) ;
            agent_info.get_end_effector_location = @(q) A.get_end_effector_location(q) ;
        end
        
        %% get collision check volume
        function out = get_collision_check_volume(A,q)
            % out_volume = A.get_collision_check_volume(q)
            %
            % Given a configuration q...
            % If the agent is 2-D, return a polyline representing the
            % agent at that configuration. This can be used with polyxpoly
            % for collision checking.
            %
            % If the agent is 3-D, return a
            % patch structure with two fields (faces and vertices) that can
            % be used with SurfaceIntersection for collision checking.
            
            if nargin < 2
                q = A.state(A.joint_state_indices,1) ;
            end
            
            [R,T] = A.get_link_rotations_and_translations(q) ;
            
            F_cell = A.collision_check_patch_data.faces ;
            V_cell = A.collision_check_patch_data.vertices ;
            
            switch A.dimension
                case 2
                    V = [] ;
                    for idx = 1:length(F_cell)
                        V_idx = R{idx}*V_cell{idx}' + T{idx} ;
                        V = [V, nan(2,1), V_idx(:,F_cell{idx})] ;
                    end
                    
                    out = V(:,2:end) ;
                case 3
                    F = [] ;
                    V = [] ;
                    
                    N_verts = 0 ;
                    
                    for idx = 1:length(F_cell)
                        F_idx = F_cell{idx} + N_verts ;
                        V_idx = (R{idx}*V_cell{idx}' + T{idx})' ;
                        
                        F = [F ; F_idx] ;
                        V = [V ; V_idx] ;
                        
                        N_verts = size(V,1) ;
                    end
                    
                    out.faces = F ;
                    out.vertices = V ;
            end
        end
        
        %% stop
        function stop(A,t)
            % A.stop(stopping_duration)
            %
            % Executes a braking maneuver from the current state for the
            % duration specified by stopping_duration; if no input is
            % given, then the arm stops for the duration A.t_stop, which is
            % 0.5 s by default.
            
            if nargin < 2
                t = A.t_stop ;
            end
            
            % time
            T_stop = [0, t] ;
            
            % desired trajectory
            Z_stop = repmat(A.state(:,end),1,2) ;
            Z_stop(A.joint_speed_indices,:) = 0 ; % joint speeds to zero
            
            % feedforward input
            joint_speed_cur = A.state(A.joint_speed_indices,end) ;
            U_brake = (zeros(size(joint_speed_cur)) - joint_speed_cur)./ t ;
            U_stop = [U_brake, zeros(size(U_brake))] ;
            
            % send command
            A.move(t,T_stop,U_stop,Z_stop) ;
        end
        
        %% forward kinematics
        function [R,T,J] = get_link_rotations_and_translations(A,time_or_config)
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
            
            if nargin < 2
                time_or_config = 0 ;
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
            
            j_locs = A.joint_locations ; % joint locations
            
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
                
                % get the value and location of the current joint
                j_idx = j_vals(idx) ;
                j_loc = j_locs(:,idx) ;
                
                % compute link rotation
                switch A.joint_types{idx}
                    case 'revolute'
                        if d == 3
                            % rotation matrix of current link
                            axis_pred = R_pred*A.joint_axes(:,idx) ;
                            R_succ = axis_angle_to_rotation_matrix_3D([axis_pred', j_idx])*R_pred ;
                        else
                            % rotation matrix of current link
                            R_succ = rotation_matrix_2D(j_idx)*R_pred ;
                        end
                        
                        % create translation
                        T_succ = T_pred + R_pred*j_loc(1:d) - R_succ*j_loc(d+1:end) ;
                    case 'prismatic'
                        % R_succ = R_pred ;
                        error('Prismatic joints are not yet supported!')
                    otherwise
                        error('Invalid joint type!')
                end
                
                % fill in rotation and translation cells
                R{s_idx} = R_succ ;
                T{s_idx} = T_succ ;
                
                % fill in the joint location (THIS IS 
                j_loc_local = j_locs((d+1):end,idx) ;
                J(:,idx) = -R_succ*j_loc_local + T_succ ;
            end
        end
        
        function  J = get_joint_locations(A,times_or_configs)
            % J = A.get_joint_locations(times_or_configs)
            %
            % Return the joint locations in 2-D or 3-D space.
            %
            % If the input is a single time t \in \R, or a single
            % configuration q \in Q, the output is a d-by-n array, where
            % n = A.n_links_and_joints and d = A.dimension.
            %
            % If the input is a 1-by-N vector of times or an n-by-N vector
            % of configurations, the output is a 1-by-N cell array where
            % each entry is the d-by-n array of joint locations for the
            % corresponding time or configuration.
            
            N = size(times_or_configs,2) ;
            if N == 1
                [~,~,J] = A.get_link_rotations_and_translations(times_or_configs) ;
            else
                J = cell(1,N) ;
                n = size(times_or_configs,1) ;
                
                if n == 1
                    % for the case of multiple times, iterate through the
                    % list and get the joint locations for each time
                    for idx = 1:N
                        [~,~,J_idx] = A.get_link_rotations_and_translations(times_or_configs(:,idx)) ;
                        J{idx} = J_idx ;
                    end
                else
                    % for the case of multiple configurations, make a cell
                    % array of the configurations and use cellfun like a
                    % heckin' matlab ninja
                    Q = mat2cell(times_or_configs, n, ones(1,N)) ;
                    J = cellfun(@(q) A.get_joint_locations_from_configuration(q),Q,'UniformOutput',false) ;
                end
            end
        end
        
        function J = get_joint_locations_from_configuration(A,q)
            % J = A.get_joint_locations_from_configuration(q)
            %
            % Return the joint locations just like A.get_joint_locations,
            % but this is faster since it does not return the rotation
            % matrices or translations of each joint, and can only take in
            % a configuration.
            
            j_vals = q ; % joint angles
            j_locs = A.joint_locations ; % joint locations
            
            % extract dimensions
            n = A.n_links_and_joints ;
            d = A.dimension ;
            
            % set up translations and rotations
            R_pred = eye(d) ;
            T_pred = zeros(d,1) ;
            
            % allocate array for the joint locations
            J = nan(d,n) ;
            
            if isa(q,'sym')
                J = sym(J) ;
            end
            
            % move through the kinematic chain and get the rotations and
            % translation of each link
            for idx = 1:n
                % get the value and location of the current joint
                j_idx = j_vals(idx) ;
                j_loc = j_locs(:,idx) ;
                
                % rotation matrix of current joint
                axis_pred = R_pred*A.joint_axes(:,idx) ;
                R_succ = axis_angle_to_rotation_matrix_3D([axis_pred', j_idx])*R_pred ;
                
                % create translation
                T_succ = T_pred + R_pred*j_loc(1:d) - R_succ*j_loc(d+1:end) ;
                
                % fill in the joint location
                j_loc_local = j_locs((d+1):end,idx) ;
                J(:,idx) = -R_succ*j_loc_local + T_succ ;
                
                % update predecessors for next iteration
                R_pred = R_succ ;
                T_pred = T_succ ;
            end
        end
        
        function J = get_end_effector_location(A,q)
            % J = A.get_end_effector_location(q)
            %
            % Does what it says! Given a configuration q, this method
            % returns the end effector location.
            
            J = A.get_joint_locations_from_configuration(q) ;
            J = J(:,end) ;
        end
        
        function [R,T,J] = forward_kinematics(A,time_or_config)
            % [R,T,J] = A.forward_kinematics(time_or_config)
            %
            % Given a time or configuration, return the link rotation
            % and translation arrays R and T, and the joint locations J.
            
            [R,T,J] = A.get_link_rotations_and_translations(time_or_config) ;
        end
        
        function ee_pos = forward_kinematics_end_effector(A,time_or_config)
            % ee_pos = forward_kinematics_end_effector(A,time_or_config)
            %
            % Return the position of the end effector as a location in
            % workspace (either 2-D or 3-D), given an input configuration;
            % the output is a vector of length A.dimension.
            [~,~,J] = A.get_link_rotations_and_translations(time_or_config) ;
            ee_pos = J(:,end) ;
        end
        
        %% inverse kinematics
        function [q,exitflag] = inverse_kinematics(A,J,q0)
            % q = A.inverse_kinematics(J)
            % q = A.inverse_kinematics(J,q0)
            % [q,exitflag] = A.inverse_kinematics(...)
            %
            % Given a desired end-effector location, or locations of all
            % the joints, J, as an A.dimension-by-(1 or A.n_links...)
            % array, attempt to find a configuration q that reaches that
            % location. This uses fmincon for nonlinear optimization of the
            % Euclidean distances squared to each joint location as the
            % cost function. The second (optional) argument is an initial
            % guess for the nonlinear solver.
            
            if nargin < 3
                q0 = A.state(A.joint_state_indices,end) ;
            end
            
            if size(J,2) > 1
                [q,exitflag] = A.inverse_kinematics_joint_locations(J,q0) ;
            else
                [q,exitflag] = A.inverse_kinematics_end_effector(J,q0) ;
            end
        end
        
        function [q,exitflag] = inverse_kinematics_joint_locations(A,J,q0)
            % q = A.inverse_kinematics_joint_locations(J)
            % q = A.inverse_kinematics_joint_locations(J,q0)
            %
            % Given joint locations J as a d-by-n array, return the
            % static configuration q as an n-by-1 vector where n is the
            % number of joints of the arm. This uses nonlinear optimization
            % (fmincon) to find the joint configuration.
            %
            % The optional second input is an initial guess for the
            % nonlinear least squares solver. If it is not returned, the
            % arm uses its last state (the column vector A.state(:,end)) as
            % the initial guess.
            
            % set up the initial config
            if nargin < 3
                q0 = A.state(A.joint_state_indices,end) ;
            end
            
            % create the least-squares function to solve for the config
            n = A.n_links_and_joints ;
            d = A.dimension ;
            opt_fun = @(x) sum(vecnorm(reshape(A.get_joint_locations(x),n*d,1) - J)) ;
            
            % create bounds on the solution
            lb = A.joint_state_limits(1,:)' ;
            ub = A.joint_state_limits(2,:)' ;
            
            % set options
            options = optimoptions('fmincon') ;
            if A.verbose < 2
                options.Display = 'off' ;
            end
            
            % optimize!
            [q,~,exitflag] = fmincon(opt_fun,q0,[],[],[],[],lb,ub,[],options) ;
        end
        
        function [q,exitflag] = inverse_kinematics_end_effector(A,J,q0)
            % q = A.inverse_kinematics_end_effector(J,q0)
            %
            % Given an end-effector location J \in \R^d where d is the
            % dimension of the agent (2 or 3), find the configuration q
            % that gets the arm to that end effector location (or fail
            % trying! aaaah!)
            
            % set up the initial guess
            if nargin < 3
                q0 = A.state(A.joint_state_indices,end) ;
            end
            
            % create the function to solve for the config
            opt_fun = @(x) sum((A.get_end_effector_location(x) - J(:)).^2) ;
            
            % create bounds on the solution
            lb = A.joint_state_limits(1,:)' ;
            ub = A.joint_state_limits(2,:)' ;
            
            % set options
            options = optimoptions('fmincon') ;
            if A.verbose < 2
                options.Display = 'off' ;
            end
            
            % optimize!
            [q,~,exitflag] = fmincon(opt_fun,q0,[],[],[],[],lb,ub,[],options) ;
        end
        
        %% moving
        function move(A,t_move,T_ref,U_ref,Z_ref)
            switch A.move_mode
                case 'integrator'
                    move@multi_link_agent(A,t_move,T_ref,U_ref,Z_ref)
                case 'direct'
                    % don't call the integrator, just assume the agent
                    % perfectly executes the reference trajectory
                    
                    % get the reference trajectory up to time t_move
                    T = 0:A.integrator_time_discretization:t_move ;
                    [U,Z] = match_trajectories(T,T_ref,U_ref,T_ref,Z_ref) ;
                    
                    % append the reference trajectory to the agent's
                    % current trajectory
                    A.commit_move_data(T,Z,T,U) ;
                otherwise
                    error('Please set A.move_mode to ''integrator'' or ''direct''!')
            end
        end
        
        function move_random(A)
            A.move(1,[0 1],zeros(A.n_inputs,2), 2.*rand(A.n_states,2) - 1)
        end
        
        %% dynamics
        function zd = dynamics(A,t,z,T,U,Z)
            % get desired torques and bound them
            u = A.LLC.get_control_inputs(A,t,z,T,U,Z) ;
            
            u = bound_array_elementwise(u,...
                A.joint_input_limits(1,:)',A.joint_input_limits(2,:)') ;
            
            % get joint speeds
            qd = z(A.joint_speed_indices) ;
            qd = bound_array_elementwise(qd,...
                A.joint_speed_limits(1,:)',A.joint_speed_limits(2,:)') ;
            
            % preallocate dynamics
            zd = zeros(A.n_states,1) ;
            
            if A.use_robotics_toolbox_model_for_dynamics_flag
                % get joint accelerations
                q = z(A.joint_state_indices) ;
                qdd = A.robotics_toolbox_model.forwardDynamics(q,qd) ;
            else
                % assume the inputs are the joint accelerations
                qdd = u(:) ;
            end
            
            % compute dynamics
            zd(A.joint_state_indices) = qd(:) ;
            zd(A.joint_speed_indices) = qdd(:) ;
        end
        
        %% integrator
        function [t_out,z_out] = integrator(A,arm_dyn,t_span,z0)
            % [tout,zout] = A.integrator(arm_dynamics,tspan,z0)
            %
            % RK4 integration with joint limits and speed limits enforced.
            
            % create time vector
            dt = A.integrator_time_discretization ;
            t_out = t_span(1):dt:t_span(end) ;
            if t_out(end) < t_span(end)
                t_out = [t_out, t_span(end)] ;
            end
            
            % preallocate trajectory output
            N_t = size(t_out,2) ;
            z_out = [z0(:), nan(A.n_states,N_t-1)] ;
            
            % run integration loop
            for tidx = 2:N_t
                % get previous state
                z_cur = z_out(:,tidx-1) ;
                t_cur = t_out(tidx-1) ;
                
                % compute RK4 terms
                k1 = arm_dyn(t_cur, z_cur) ;
                k2 = arm_dyn(t_cur + dt/2, z_cur + dt*k1/2) ;
                k3 = arm_dyn(t_cur + dt/2, z_cur + dt*k2/2) ;
                k4 = arm_dyn(t_cur + dt, z_cur + dt*k3) ;
                
                % compute summed term
                dzdt = (1/6)*(k1 + 2*k2 + 2*k3 + k4) ;
                
                % compute next state
                z_new = z_cur + dt*dzdt ;
                
                % apply state limits
                joint_values = z_new(A.joint_state_indices)' ;
                joint_values = max([joint_values ; A.joint_state_limits(1,:)],[],1) ;
                joint_values = min([joint_values ; A.joint_state_limits(2,:)],[],1) ;
                z_new(A.joint_state_indices) = joint_values ;
                
                % apply speed limits
                joint_speeds = z_new(A.joint_speed_indices)' ;
                joint_speeds = max([joint_speeds; A.joint_speed_limits(1,:)],[],1) ;
                joint_speeds = min([joint_speeds; A.joint_speed_limits(2,:)],[],1) ;
                z_new(A.joint_speed_indices) = joint_speeds ;
                
                % save new state
                z_out(:,tidx) = z_new(:) ;
            end
        end
        
        %% plotting
        function plot(A,~)
            A.plot_at_time(A.time(end)) ;
        end
        
        function plot_at_time(A,t,~)
            if nargin < 2
                t = 0 ;
            end
            
            A.plot_baselink() ;
            A.plot_links(t) ;
        end
        
        function plot_baselink(A)
            % plot baselink
            BF = A.link_plot_data.baselink_faces ;
            BV = A.link_plot_data.baselink_vertices ;
            
            if check_if_plot_is_available(A,'baselink')
                A.plot_data.baselink.Faces = BF ;
                A.plot_data.baselink.Vertices = BV ;
            else
                baselink_data = patch('Faces',BF,'Vertices',BV,...
                    'FaceColor',A.link_plot_face_color,...
                    'FaceAlpha',A.link_plot_face_opacity,...
                    'EdgeColor',A.link_plot_edge_color,...
                    'LineWidth',A.link_plot_edge_width,...
                    'EdgeAlpha',A.link_plot_edge_opacity) ;
                A.plot_data.baselink = baselink_data ;
            end
        end
        
        function plot_links(A,time_or_config)
            % get the rotations and translations at the current time
            [R,T] = A.get_link_rotations_and_translations(time_or_config) ;
            
            % generate plot data for each link
            link_verts = cell(1,A.n_links_and_joints) ;
            for idx = 1:A.n_links_and_joints
                link_verts{idx} = (R{idx}*A.link_plot_data.link_vertices{idx}' + ...
                    T{idx})' ;
            end
            
            if check_if_plot_is_available(A,'links')
                for idx = 1:A.n_links_and_joints
                    A.plot_data.links(idx).Faces = A.link_plot_data.link_faces{idx} ;
                    A.plot_data.links(idx).Vertices = link_verts{idx} ;
                end
            else
                link_array = [] ;
                for idx = 1:A.n_links_and_joints
                    link_data = patch('Faces',A.link_plot_data.link_faces{idx},...
                        'Vertices',link_verts{idx},...
                        'FaceColor',A.link_plot_face_color,...
                        'FaceAlpha',A.link_plot_face_opacity,...
                        'EdgeColor',A.link_plot_edge_color,...
                        'LineWidth',A.link_plot_edge_width,...
                        'EdgeAlpha',A.link_plot_edge_opacity) ;
                    link_array = [link_array, link_data] ;
                end
                A.plot_data.links = link_array ;
            end
        end
        
        function lims = get_axis_lims(A)
            % figure out the maximum length of the arm
            L = sum(A.link_sizes(1,:)) ;
            
            % create axis limits
            switch A.dimension
                case 2
                    lims = [-L,L,-L,L] ;
                case 3
%                     switch A.floor_normal_axis
%                         case 1
%                             lims = [-L, L, -L, L, -L, L] ;
%                         case 2
%                             lims = [-L, L, 0, L, -L, L] ;
%                         case 3
%                             lims = [-L,L,-L,L,0,L] ;
%                     end
                    lims = [-L, L, -L, L, -L, L] ;
                    
                    % in case base of robot is not at [0;0;0]:
                    lims = lims + [A.joint_locations(1, 1)*ones(1, 2),...
                        A.joint_locations(2, 1)*ones(1, 2),...
                        A.joint_locations(3, 1)*ones(1, 2)];
            end
        end
    end
end