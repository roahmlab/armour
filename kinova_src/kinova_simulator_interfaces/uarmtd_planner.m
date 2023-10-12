classdef uarmtd_planner < robot_arm_generic_planner
    % UARMTD_PLANNER implements new polynomial zonotope collision-avoidance
    % and input constraints
    
    properties
        % housekeeping
        time_discretization = 0.01;        
        trajopt_start_tic;
        iter = 0;
        first_iter_pause_flag = true ;
        
        % jrs and agent info
        jrs_info;
        agent_info;

        % cost type
        use_q_plan_for_cost = false; % otherwise use q_stop (q at final time)

        % constraints:
        constraints = {}; % cell will contain functions of $k$ for evaluating constraints
        grad_constraints = {}; % cell will contain functions of $k$ for evaluating gradients
        obs_constraints = {};
        grad_obs_constraints = {};

        % smooth obstacle constraints:
        smooth_obs_constraints;
        smooth_obs_constraints_A;
        smooth_obs_lambda_index;
        
        % for turning on/off constraint types
        save_FO_zono_flag = true;
        input_constraints_flag = true;
        grasp_constraints_flag = false; % must have input constraints turned on!
        smooth_obstacle_constraints_flag = false;
        use_robust_input = true; % turn this off to only consider nominal passivity-based controller

        % for obstacle avoidance:
        combs;

        % for JRSs and trajectories:
        taylor_degree = 1;
        traj_type = 'bernstein'; % choose 'orig' for original ARMTD, or 'bernstein'
%         use_waypoint_for_bernstein_center = false; % if true, centers bernstein final configuration range around waypoint
        jrs_type = 'online';
        use_cuda = false;

        kinova_test_folder_path = '';
    end
    
    methods
        function P = uarmtd_planner(varargin)
            t_move = 0.5;
            HLP = robot_arm_straight_line_HLP( );
            P@robot_arm_generic_planner('t_move', t_move, 'HLP', HLP, ...
                varargin{:}) ;
            
            % hard code planning time...
            P.t_plan = 0.5;

            % init info object
            P.init_info()
            
            % initialize combinations (for obstacle avoidance constraints)
            P.combs.maxcombs = 200;
            P.combs.combs = generate_combinations_upto(200);

            data = load('c_kvi.mat');
            P.jrs_info.c_kvi = data.c_kvi;

            data = load('kinova_test_folder_path');
            P.kinova_test_folder_path = data.kinova_test_folder_path;
        end
        
        function init_info(P)
            P.info = struct('T',[],'U',[],'Z',[],'waypoint',[],...
                'obstacles',[],'q_0',[],'q_dot_0',[],'k_opt',[],...
                'desired_trajectory', [], 't_move', [], ...
                'FO_zono', [], 'sliced_FO_zono', []) ;
        end
        
        function [T, U, Z, info] = replan(P,agent_info,world_info)
            P.vdisp('Replanning!',5)
            
            % get current state of robot
            P.agent_info = agent_info;

            q_0 = agent_info.reference_state(P.arm_joint_state_indices, end) ;
            q_dot_0 = agent_info.reference_state(P.arm_joint_speed_indices, end) ;
            q_ddot_0 = agent_info.reference_acceleration(:, end);
            % q_ddot_0 = zeros(size(q_0)); % need to pass this in for bernstein!!!

            if ~P.use_cuda
                % if bounds are +- inf, set to small/large number
                joint_limit_infs = isinf(agent_info.joint_state_limits) ;
                speed_limit_infs = isinf(agent_info.joint_speed_limits) ;
                input_limit_infs = isinf(agent_info.joint_input_limits) ;
    
                P.agent_info.joint_state_limits(1,joint_limit_infs(1,:)) = -200*pi ;
                P.agent_info.joint_state_limits(2,joint_limit_infs(2,:)) = +200*pi ;            
                P.agent_info.joint_speed_limits(1,speed_limit_infs(1,:)) = -200*pi ;
                P.agent_info.joint_speed_limits(2,speed_limit_infs(2,:)) = +200*pi ;
                P.agent_info.joint_input_limits(1,input_limit_infs(1,:)) = -200*pi ;
                P.agent_info.joint_input_limits(2,input_limit_infs(2,:)) = +200*pi ;
    
                % generate a waypoint in configuration space
                P.vdisp('Generating cost function',6)
                if P.first_iter_pause_flag && P.iter == 0
                   pause; 
                end
                P.iter = P.iter + 1;
                planning_time = tic;
                q_des = P.HLP.get_waypoint(agent_info,world_info,P.lookahead_distance) ;
                if isempty(q_des)
                    P.vdisp('Waypoint creation failed! Using global goal instead.', 3)
                    q_des = P.HLP.goal ;
                end
                            
                % get current obstacles and create constraints
                P.vdisp('Generating constraints',6)
                [P, FO] = P.generate_constraints(q_0, q_dot_0, q_ddot_0, world_info.obstacles);
                
                % optimize
                P.vdisp('Replan is calling trajopt!',8)
                [k_opt, trajopt_failed] = P.trajopt(q_0, q_dot_0, q_ddot_0, q_des);
                if P.smooth_obstacle_constraints_flag
                    k_opt = k_opt(1:P.agent_info.params.pz_nominal.num_q);
                end
    
                % process result            
                if ~trajopt_failed
                    P.vdisp('New trajectory found!',3);
                else % no safe trajectory parameter found:
                    P.vdisp('Unable to find new trajectory!',3)
                    k_opt = nan;
                end
                toc(planning_time);
            else
                if strcmp(P.traj_type, 'bernstein')
                    P.jrs_info.n_t = 96;
                    P.jrs_info.n_q = 7;
                    P.jrs_info.n_k = 7;
                    P.jrs_info.c_k_bernstein = zeros(7,1);
                    % !!!!!!
                    % Make sure this is consistent with the k_range in
                    % kinova_src/kinova_simulator_interfaces/kinova_planner_realtime/Parameters.h 
                    % !!!!!!
%                     P.jrs_info.g_k_bernstein = [pi/24; pi/24; pi/24; pi/24; pi/24; pi/24; pi/24];
                    P.jrs_info.g_k_bernstein = [pi/48; pi/48; pi/48; pi/48; pi/48; pi/48; pi/48];
    
                    q_des = P.HLP.get_waypoint(agent_info,world_info,P.lookahead_distance) ;
                    if isempty(q_des)
                        P.vdisp('Waypoint creation failed! Using global goal instead.', 3)
                        q_des = P.HLP.goal ;
                    end
    
                    % organize input to cuda program
                    P.vdisp('Calling CUDA & C++ Program!',3);
                    cuda_input_file = fopen([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime/buffer/armour.in'], 'w');
                
                    for ind = 1:length(q_0)
                        fprintf(cuda_input_file, '%.10f ', q_0(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    for ind = 1:length(q_dot_0)
                        fprintf(cuda_input_file, '%.10f ', q_dot_0(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    for ind = 1:length(q_ddot_0)
                        fprintf(cuda_input_file, '%.10f ', q_ddot_0(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    for ind = 1:length(q_des)
                        fprintf(cuda_input_file, '%.10f ', q_des(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    fprintf(cuda_input_file, '%d\n', max(length(world_info.obstacles), 0));
                    for obs_ind = 1:length(world_info.obstacles)
                        temp = reshape(world_info.obstacles{obs_ind}.Z, [1,size(world_info.obstacles{obs_ind}.Z,1) * size(world_info.obstacles{obs_ind}.Z,2)]);
                        for ind = 1:length(temp)
                            fprintf(cuda_input_file, '%.10f ', temp(ind));
                        end
                        fprintf(cuda_input_file, '\n');
                    end
    
                    fclose(cuda_input_file);
    
                    % call cuda program in terminal
                    % you have to be in the proper path!
                    terminal_output = system('./../kinova_simulator_interfaces/kinova_planner_realtime/armour_main');
    
                    if terminal_output == 0
                        data = readmatrix('armour.out', 'FileType', 'text');
                        k_opt = data(1:end-1);
                        planning_time = data(end) / 1000.0; % original data is milliseconds
    
                        if length(k_opt) == 1
                            P.vdisp('Unable to find new trajectory!',3)
                            k_opt = nan;
                        else
                            P.vdisp('New trajectory found!',3);
                            for i = 1:length(k_opt)
                                fprintf('%7.6f ', k_opt(i));
                            end
                            fprintf('\n');
                        end
                    else
                        error('CUDA program error! Check the executable path in armour-dev/kinova_src/kinova_simulator_interfaces/uarmtd_planner');
                    end
    
                    if terminal_output == 0
                        % read FRS information if needed
                        link_frs_center = readmatrix('armour_joint_position_center.out', 'FileType', 'text');
                        link_frs_generators = readmatrix('armour_joint_position_radius.out', 'FileType', 'text');
                        if P.input_constraints_flag
                            control_input_radius = readmatrix('armour_control_input_radius.out', 'FileType', 'text');
                        else
                            control_input_radius = [];
                        end
                        constraints_value = readmatrix('armour_constraints.out', 'FileType', 'text');

                        record_tids = [1:10:P.jrs_info.n_t, P.jrs_info.n_t];
                        link_frs_vertices = cell(7,length(record_tids));
                        for tid = record_tids
                            for j = 1:agent_info.n_links_and_joints
                                c = link_frs_center((tid-1)*7+j, :)';
                                g = link_frs_generators( ((tid-1)*7+j-1)*3+1 : ((tid-1)*7+j)*3, :);
                                Z = zonotope(c, g);
                                link_frs_vertices{j,tid} = [link_frs_vertices{j}; vertices(Z)'];
                            end
                        end
                    else
                        k_opt = nan;
                    end
                elseif strcmp(P.traj_type, 'orig')
                    P.jrs_info.n_t = 100;
                    P.jrs_info.n_q = 7;
                    P.jrs_info.n_k = 7;
                    P.jrs_info.c_k_orig = zeros(7,1);
                    P.jrs_info.g_k_orig = zeros(7,1);

                    q_des = P.HLP.get_waypoint(agent_info,world_info,P.lookahead_distance) ;
                    if isempty(q_des)
                        P.vdisp('Waypoint creation failed! Using global goal instead.', 3)
                        q_des = P.HLP.goal ;
                    end

                    % organize input to cuda program
                    P.vdisp('Calling CUDA & C++ Program!',3);
                    JRS_filename = cell(P.jrs_info.n_q,1);
                    for i = 1:length(q_0)
                        [~, closest_idx] = min(abs(q_dot_0(i) - P.jrs_info.c_kvi));
                        filename = sprintf('%sJRS_%0.3f.mat', 'offline_jrs/orig_parameterization/', P.jrs_info.c_kvi(closest_idx));
                        
                        JRS_filename{i} = filename;
                    end
    
                    cuda_input_file = fopen([P.kinova_test_folder_path, '/kinova_simulator_interfaces/kinova_planner_realtime_armtd_comparison/buffer/armtd.in'], 'w');
                    
                    % q_0
                    for ind = 1:length(q_0)
                        fprintf(cuda_input_file, '%.10f ', q_0(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    
                    % q_dot_0
                    for ind = 1:length(q_dot_0)
                        fprintf(cuda_input_file, '%.10f ', q_dot_0(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    
                    % q_des
                    for ind = 1:length(q_des)
                        fprintf(cuda_input_file, '%.10f ', q_des(ind));
                    end
                    fprintf(cuda_input_file, '\n');
                    
                    for ind = 1:P.jrs_info.n_k
                        data = load(JRS_filename{ind});
                    
                        % cos(q)
                        for i = 1:P.jrs_info.n_t
                            c = data.JRS{i}.Z(1,1);
                            fprintf(cuda_input_file, '%.10f ', c);
                        end
                        fprintf(cuda_input_file, '\n');
                        for i = 1:P.jrs_info.n_t
                            g = data.JRS{i}.Z(1,2);
                            fprintf(cuda_input_file, '%.10f ', g);
                        end
                        fprintf(cuda_input_file, '\n');
                        for i = 1:P.jrs_info.n_t
                            r = sum(abs(data.JRS{i}.Z(1,3:end)), 2);
                            fprintf(cuda_input_file, '%.10f ', r);
                        end
                        fprintf(cuda_input_file, '\n');
                        
                        % sin(q)
                        for i = 1:P.jrs_info.n_t
                            c = data.JRS{i}.Z(2,1);
                            fprintf(cuda_input_file, '%.10f ', c);
                        end
                        fprintf(cuda_input_file, '\n');
                        for i = 1:P.jrs_info.n_t
                            g = data.JRS{i}.Z(2,2);
                            fprintf(cuda_input_file, '%.10f ', g);
                        end
                        fprintf(cuda_input_file, '\n');
                        for i = 1:P.jrs_info.n_t
                            r = sum(abs(data.JRS{i}.Z(2,3:end)), 2);
                            fprintf(cuda_input_file, '%.10f ', r);
                        end
                        fprintf(cuda_input_file, '\n');
                        
                        P.jrs_info.g_k_orig(ind) = data.JRS{1}.Z(5,2);
                        fprintf(cuda_input_file, '%.10f\n ', P.jrs_info.g_k_orig(ind));
                    end

                    fprintf(cuda_input_file, '%d\n', max(length(world_info.obstacles), 0));
                    for obs_ind = 1:length(world_info.obstacles)
                        temp = reshape(world_info.obstacles{obs_ind}.Z, [1,size(world_info.obstacles{obs_ind}.Z,1) * size(world_info.obstacles{obs_ind}.Z,2)]);
                        for ind = 1:length(temp)
                            fprintf(cuda_input_file, '%.10f ', temp(ind));
                        end
                        fprintf(cuda_input_file, '\n');
                    end
    
                    fclose(cuda_input_file);

                    % call cuda program in terminal
                    % you have to be in the proper path!
                    terminal_output = system('./../kinova_simulator_interfaces/kinova_planner_realtime_armtd_comparison/armtd_main');
    
                    if terminal_output == 0
                        data = readmatrix('armtd.out', 'FileType', 'text');
                        k_opt = data(1:end-1);
                        planning_time = data(end) / 1000.0; % original data is milliseconds
    
                        if length(k_opt) == 1
                            P.vdisp('Unable to find new trajectory!',3)
                            k_opt = nan;
                        else
                            P.vdisp('New trajectory found!',3);
                            for i = 1:length(k_opt)
                                fprintf('%7.6f ', k_opt(i));
                            end
                            fprintf('\n');
                        end
                    else
                        error('CUDA program error! Check the executable path in armour-dev/kinova_src/kinova_simulator_interfaces/uarmtd_planner');
                    end
    
                    if terminal_output == 0
                        % read FRS information if needed
                        link_frs_center = readmatrix('armtd_joint_position_center.out', 'FileType', 'text');
                        link_frs_generators = readmatrix('armtd_joint_position_radius.out', 'FileType', 'text');
                        constraints_value = readmatrix('armtd_constraints.out', 'FileType', 'text');

                        record_tids = [1:10:P.jrs_info.n_t, P.jrs_info.n_t];
                        link_frs_vertices = cell(agent_info.n_links_and_joints,length(record_tids));
                        for i = 1:length(record_tids)
                            tid = record_tids(i);
                            for j = 1:agent_info.n_links_and_joints
                                c = link_frs_center((tid-1)*agent_info.n_links_and_joints+j, :)';
                                g = link_frs_generators( ((tid-1)*agent_info.n_links_and_joints+j-1)*3+1 : ((tid-1)*agent_info.n_links_and_joints+j)*3, :);
                                Z = zonotope(c, g);
                                link_frs_vertices{j,i} = vertices(Z)';
                            end
                        end
                    else
                        k_opt = nan;
                    end
                else
                    error('Unrecognized trajectory type!');
                end
            end

            disturbance = 2 * rand - 1;
            disturbance = 0.02 * disturbance / norm(disturbance);
            q_0 = q_0 + disturbance;
            disturbance = 2 * rand - 1;
            disturbance = 0.02 * disturbance / norm(disturbance);
            q_dot_0 = q_dot_0 + disturbance;
            
            % save info
            P.info.planning_time = planning_time;
            P.info.desired_trajectory = [P.info.desired_trajectory, {@(t) P.desired_trajectory(q_0, q_dot_0, q_ddot_0, t, k_opt)}];
            P.info.t_move = [P.info.t_move, {P.t_move}];
            P.info.waypoint = [P.info.waypoint, {q_des}] ;
            P.info.obstacles = [P.info.obstacles, {world_info.obstacles}] ;
            P.info.q_0 = [P.info.q_0, {q_0}] ;
            P.info.q_dot_0 = [P.info.q_dot_0, {q_dot_0}] ;
            P.info.k_opt = [P.info.k_opt, {k_opt}] ;
            if P.save_FO_zono_flag
                if ~P.use_cuda
                    for i = 1:P.jrs_info.n_t
                        for j = 1:P.agent_info.params.pz_nominal.num_bodies
                            FO_zono{i}{j} = zonotope(FO{i}{j});
                            if trajopt_failed
                                % no safe slice
                                sliced_FO_zono{i}{j} = [];
                            else
                                % slice and save
                                fully_sliceable_tmp = polyZonotope_ROAHM(FO{i}{j}.c, FO{i}{j}.G, [], FO{i}{j}.expMat, FO{i}{j}.id);
                                sliced_FO_zono{i}{j} = zonotope([slice(fully_sliceable_tmp, k_opt), FO{i}{j}.Grest]);
                            end
                        end
                    end
                    P.info.FO_zono = [P.info.FO_zono, {FO_zono}];
                    P.info.sliced_FO_zono = [P.info.sliced_FO_zono, {sliced_FO_zono}];
                else
                    P.info.sliced_FO_zono = [P.info.sliced_FO_zono, {link_frs_vertices}]; % disable recording for now
                end
            end

            % create outputs:
            T = 0:P.time_discretization:P.t_stop ;
            U = zeros(agent_info.n_inputs, length(T));
            Z = zeros(agent_info.n_states, length(T));
            for i = 1:length(T)
                if ~isempty(P.info.desired_trajectory{end})
                    [q_tmp, qd_tmp, ~] = P.info.desired_trajectory{end}(T(i));
                else
                    q_tmp = nan;
                    qd_tmp = nan;
                end
                Z(agent_info.joint_state_indices, i) = q_tmp;
                Z(agent_info.joint_speed_indices, i) = qd_tmp;
            end
            info = P.info;

        end
        
        function [P, FO] = generate_constraints(P, q_0, q_dot_0, q_ddot_0, O, waypoint)
            %% get JRSs:
            if strcmp(P.jrs_type, 'offline')
                if ~strcmp(P.traj_type, 'orig')
                    error('Offline JRS only implemented for original parameterization');
                end
                [q_des, dq_des, ddq_des, q, dq, dq_a, ddq_a, R_des, R_t_des, R, R_t, jrs_info] = load_offline_jrs(q_0, q_dot_0, q_ddot_0,...
                    P.agent_info.params.pz_nominal.joint_axes, P.use_robust_input);
            else
%                 if P.use_waypoint_for_bernstein_center
%                     [q_des, dq_des, ddq_des, q, dq, dq_a, ddq_a, R_des, R_t_des, R, R_t, jrs_info] = create_jrs_online(q_0, q_dot_0, q_ddot_0,...
%                         P.agent_info.params.pz_nominal.joint_axes, P.taylor_degree, P.traj_type, P.use_robust_input, waypoint);
%                 else
                    [q_des, dq_des, ddq_des, q, dq, dq_a, ddq_a, R_des, R_t_des, R, R_t, jrs_info] = create_jrs_online(q_0, q_dot_0, q_ddot_0,...
                        P.agent_info.params.pz_nominal.joint_axes, P.taylor_degree, P.traj_type, P.use_robust_input, P.agent_info.LLC_info);
%                 end
            end
            P.jrs_info = jrs_info;

            %% create FO and input poly zonotopes:
            % set up zeros and overapproximation of r
            for j = 1:jrs_info.n_q
                zero_cell{j, 1} = polyZonotope_ROAHM(0);
                if P.use_robust_input
                    r{j, 1} = polyZonotope_ROAHM(0, [], P.agent_info.LLC_info.ultimate_bound);
                end
            end
            
            % get forward kinematics and forward occupancy
            for i = 1:jrs_info.n_t
               [R_w{i, 1}, p_w{i, 1}] = pzfk(R{i, 1}, P.agent_info.params.pz_nominal); 
               for j = 1:P.agent_info.params.pz_nominal.num_bodies
                  FO{i, 1}{j, 1} = R_w{i, 1}{j, 1}*P.agent_info.link_poly_zonotopes{j, 1} + p_w{i, 1}{j, 1}; 
                  FO{i, 1}{j, 1} = reduce(FO{i, 1}{j, 1}, 'girard', P.agent_info.params.pz_interval.zono_order);
                  FO{i, 1}{j, 1} = remove_dependence(FO{i, 1}{j, 1}, jrs_info.k_id(end));
               end
            end

            % get nominal inputs, disturbances, possible lyapunov functions:
            % can possibly speed up by using tau_int - tau_int to get disturbance.
            if P.input_constraints_flag
                for i = 1:jrs_info.n_t
                    tau_nom{i, 1} = poly_zonotope_rnea(R{i}, R_t{i}, dq{i}, dq_a{i}, ddq_a{i}, true, P.agent_info.params.pz_nominal);
                    if P.use_robust_input
                        [tau_int{i, 1}, f_int{i, 1}, n_int{i, 1}] = poly_zonotope_rnea(R{i}, R_t{i}, dq{i}, dq_a{i}, ddq_a{i}, true, P.agent_info.params.pz_interval);
                        for j = 1:jrs_info.n_q
                            w{i, 1}{j, 1} = tau_int{i, 1}{j, 1} - tau_nom{i, 1}{j, 1};
                            w{i, 1}{j, 1} = reduce(w{i, 1}{j, 1}, 'girard', P.agent_info.params.pz_interval.zono_order);
                        end
                        V_cell = poly_zonotope_rnea(R{i}, R_t{i}, zero_cell, zero_cell, r, false, P.agent_info.params.pz_interval);
                        V{i, 1} = 0;
                        for j = 1:jrs_info.n_q
                            V{i, 1} = V{i, 1} + 0.5.*r{j, 1}.*V_cell{j, 1};
                            V{i, 1} = reduce(V{i, 1}, 'girard', P.agent_info.params.pz_interval.zono_order);
                        end
                        V_diff{i, 1} = V{i, 1} - V{i, 1};
                        V_diff{i, 1} = reduce(V_diff{i, 1}, 'girard', P.agent_info.params.pz_interval.zono_order);
                        V_diff_int{i, 1} = interval(V_diff{i, 1});
                    end
                end

                % % get max effect of disturbance r'*w... couple of ways to do this
                % % can overapproximate r and compute r'*w as a PZ
                % % can slice w first, then get the norm?
                % for i = 1:jrs_info.n_t
                %     r_dot_w{i, 1} = 0;
                %     for j = 1:jrs_info.n_q
                %         r_dot_w{i, 1} = r_dot_w{i, 1} + r{j, 1}.*w{i, 1}{j, 1};
                %         r_dot_w{i, 1} = reduce(r_dot_w{i, 1}, 'girard', P.params.pz_interval.zono_order);
                %     end
                % end
                % can get ||w|| <= ||\rho(\Phi)||, and compute the norm using interval arithmetic
                if P.use_robust_input
                    for i = 1:jrs_info.n_t
                        for j = 1:jrs_info.n_q
                            w_int{i, 1}(j, 1) = interval(w{i, 1}{j, 1});
                        end
                        rho_max{i, 1} = norm(max(abs(w_int{i, 1}.inf), abs(w_int{i, 1}.sup)));
                    end

                    % compute robust input bound tortatotope:
                    for i = 1:jrs_info.n_t
                        v_norm{i, 1} = (P.agent_info.LLC_info.alpha_constant*V_diff_int{i, 1}.sup).*(1/P.agent_info.LLC_info.ultimate_bound) + rho_max{i, 1};
    %                     v_norm{i, 1} = reduce(v_norm{i, 1}, 'girard', P.agent_info.params.pz_interval.zono_order);
                    end
                else
                    for i = 1:jrs_info.n_t
                        v_norm{i, 1} = 0;
                    end
                end

                % compute total input tortatotope
                for i = 1:jrs_info.n_t
                    for j = 1:jrs_info.n_q
                        u_ub_tmp = tau_nom{i, 1}{j, 1} + v_norm{i, 1};
                        u_lb_tmp = tau_nom{i, 1}{j, 1} - v_norm{i, 1};
                        u_ub_tmp = remove_dependence(u_ub_tmp, jrs_info.k_id(end));
                        u_lb_tmp = remove_dependence(u_lb_tmp, jrs_info.k_id(end));
                        u_ub_buff = sum(abs(u_ub_tmp.Grest));
                        u_lb_buff = -sum(abs(u_lb_tmp.Grest));
                        u_ub{i, 1}{j, 1} = polyZonotope_ROAHM(u_ub_tmp.c + u_ub_buff, u_ub_tmp.G, [], u_ub_tmp.expMat, u_ub_tmp.id) - P.agent_info.joint_input_limits(2, j);
                        u_lb{i, 1}{j, 1} = -1*polyZonotope_ROAHM(u_lb_tmp.c + u_lb_buff, u_lb_tmp.G, [], u_lb_tmp.expMat, u_lb_tmp.id) + P.agent_info.joint_input_limits(1, j);
                    end
                end
                
                % grasp constraints:
                if P.grasp_constraints_flag
                    % ASSUMING SURFACE NORMAL IS POSITIVE Z-DIRECTION
                    % fill this out still
                end
            end
            
            % joint limit constraint setup
            for i = 1:jrs_info.n_t
                for j = 1:jrs_info.n_q
                    q_lim_tmp = q{i, 1}{j, 1};
                    dq_lim_tmp = dq{i, 1}{j, 1};
                    q_lim_tmp = remove_dependence(q_lim_tmp, jrs_info.k_id(end));
                    dq_lim_tmp = remove_dependence(dq_lim_tmp, jrs_info.k_id(end));
                    q_buf = sum(abs(q_lim_tmp.Grest));
                    dq_buf = sum(abs(dq_lim_tmp.Grest));
                    q_ub{i, 1}{j, 1} = polyZonotope_ROAHM(q_lim_tmp.c + q_buf, q_lim_tmp.G, [], q_lim_tmp.expMat, q_lim_tmp.id) - P.agent_info.joint_state_limits(2, j);
                    q_lb{i, 1}{j, 1} = -1*polyZonotope_ROAHM(q_lim_tmp.c + q_buf, q_lim_tmp.G, [], q_lim_tmp.expMat, q_lim_tmp.id) + P.agent_info.joint_state_limits(1, j);
                    dq_ub{i, 1}{j, 1} = polyZonotope_ROAHM(dq_lim_tmp.c + dq_buf, dq_lim_tmp.G, [], dq_lim_tmp.expMat, dq_lim_tmp.id) - P.agent_info.joint_speed_limits(2, j);
                    dq_lb{i, 1}{j, 1} = -1*polyZonotope_ROAHM(dq_lim_tmp.c + dq_buf, dq_lim_tmp.G, [], dq_lim_tmp.expMat, dq_lim_tmp.id) + P.agent_info.joint_speed_limits(1, j);
                end
            end

            %% start making constraints
            P.constraints = {};
            P.grad_constraints = {};
            if ~P.smooth_obstacle_constraints_flag
                P.obs_constraints = {};
                P.grad_obs_constraints = {};
            else
                P.smooth_obs_constraints = {};
                P.smooth_obs_constraints_A = {};
                P.smooth_obs_lambda_index = {};
            end

            % obstacle avoidance constraints
            for i = 1:jrs_info.n_t
                for j = 1:P.agent_info.params.pz_nominal.num_bodies
                    for o = 1:length(O) % for each obstacle
                        
                        % first, check if constraint is necessary
                        O_buf = [O{o}.Z, FO{i, 1}{j, 1}.G, FO{i, 1}{j, 1}.Grest];
                        [A_obs, b_obs] =  polytope_PH(O_buf, P.combs); % get polytope form
                        if ~(all(A_obs*FO{i, 1}{j, 1}.c - b_obs <= 0, 1))
                            continue;
                        end
                        
                        % reduce FO so that polytope_PH has fewer
                        % directions to consider
                        FO{i, 1}{j, 1} = reduce(FO{i, 1}{j, 1}, 'girard', 3);
                        
                        % now create constraint
                        FO_buf = FO{i, 1}{j, 1}.Grest; % will buffer by non-sliceable gens
                        O_buf = [O{o}.Z, FO_buf]; % describes buffered obstacle zonotope
                        [A_obs, b_obs] = polytope_PH(O_buf, P.combs); % get polytope form

                        % constraint PZ:
                        FO_tmp = polyZonotope_ROAHM(FO{i, 1}{j, 1}.c, FO{i, 1}{j, 1}.G, [], FO{i, 1}{j, 1}.expMat, FO{i, 1}{j, 1}.id);
                        obs_constraint_pz = A_obs*FO_tmp - b_obs;

                        % turn into function
                        obs_constraint_pz_slice = @(k) slice(obs_constraint_pz, k);

                        % add gradients
                        grad_obs_constraint_pz = grad(obs_constraint_pz, P.jrs_info.n_q);
                        grad_obs_constraint_pz_slice = @(k) cellfun(@(C) slice(C, k), grad_obs_constraint_pz, 'UniformOutput', false);
                        
                        % save
                        if ~P.smooth_obstacle_constraints_flag
                            P.obs_constraints{end+1, 1} = @(k) P.evaluate_obs_constraint(obs_constraint_pz_slice, grad_obs_constraint_pz_slice, k);
                        else
                            P.smooth_obs_constraints_A{end+1, 1} = A_obs;
                            P.smooth_obs_constraints{end+1, 1} = @(k, lambda) P.evaluate_smooth_obs_constraint(obs_constraint_pz_slice, grad_obs_constraint_pz_slice, k, lambda);
                            if isempty(P.smooth_obs_lambda_index)
                                P.smooth_obs_lambda_index{end+1, 1} = (1:size(obs_constraint_pz.c, 1))';
                            else
                                P.smooth_obs_lambda_index{end+1, 1} = P.smooth_obs_lambda_index{end}(end, 1) + (1:size(obs_constraint_pz.c, 1))';
                            end
                        end
                    end
                end
            end 

            % input constraints
            if P.input_constraints_flag
                for i = 1:jrs_info.n_t
                    for j = 1:jrs_info.n_q

                        % first check if constraints are necessary, then add
                        u_ub_int = interval(u_ub{i, 1}{j, 1});
                        if ~(u_ub_int.sup < 0)
                            % no constraint necessary
%                             continue;
%                         else
                            % add constraint and gradient
                            fprintf('ADDED UPPER BOUND INPUT CONSTRAINT ON JOINT %d AT TIME %d \n', j, i);
                            P.constraints{end+1, 1} = @(k) slice(u_ub{i, 1}{j, 1}, k);
                            grad_u_ub = grad(u_ub{i, 1}{j, 1}, P.jrs_info.n_q);
                            P.grad_constraints{end+1, 1} = @(k) cellfun(@(C) slice(C, k), grad_u_ub);
                        end

                        u_lb_int = interval(u_lb{i, 1}{j, 1});
                        if ~(u_lb_int.sup < 0)
                            % no constraint necessary
%                             continue;
%                         else
                            % add constraint and gradient
                            fprintf('ADDED LOWER BOUND INPUT CONSTRAINT ON JOINT %d AT TIME %d \n', j, i);
                            P.constraints{end+1, 1} = @(k) slice(u_lb{i, 1}{j, 1}, k);
                            grad_u_lb = grad(u_lb{i, 1}{j, 1}, P.jrs_info.n_q);
                            P.grad_constraints{end+1, 1} = @(k) cellfun(@(C) slice(C, k), grad_u_lb);
                        end
                    end
                end
            end
            
            % joint limit constraints
            for i = 1:jrs_info.n_t
                for j = 1:jrs_info.n_q
                    % check if constraint necessary, then add
                    q_ub_int = interval(q_ub{i, 1}{j, 1});
                    if ~(q_ub_int.sup < 0)
                        fprintf('ADDED UPPER BOUND JOINT POSITION CONSTRAINT ON JOINT %d AT TIME %d \n', j, i);
                        P.constraints{end+1, 1} = @(k) slice(q_ub{i, 1}{j, 1}, k);
                        grad_q_ub = grad(q_ub{i, 1}{j, 1}, P.jrs_info.n_q);
                        P.grad_constraints{end+1, 1} = @(k) cellfun(@(C) slice(C, k), grad_q_ub);
                    end
                    
                    q_lb_int = interval(q_lb{i, 1}{j, 1});
                    if ~(q_lb_int.sup < 0)
                        fprintf('ADDED LOWER BOUND JOINT POSITION CONSTRAINT ON JOINT %d AT TIME %d \n', j, i);
                        P.constraints{end+1, 1} = @(k) slice(q_lb{i, 1}{j, 1}, k);
                        grad_q_lb = grad(q_lb{i, 1}{j, 1}, P.jrs_info.n_q);
                        P.grad_constraints{end+1, 1} = @(k) cellfun(@(C) slice(C, k), grad_q_lb);
                    end
                    
                    dq_ub_int = interval(dq_ub{i, 1}{j, 1});
                    if ~(dq_ub_int.sup < 0)
                        fprintf('ADDED UPPER BOUND JOINT VELOCITY CONSTRAINT ON JOINT %d AT TIME %d \n', j, i);
                        P.constraints{end+1, 1} = @(k) slice(dq_ub{i, 1}{j, 1}, k);
                        grad_dq_ub = grad(dq_ub{i, 1}{j, 1}, P.jrs_info.n_q);
                        P.grad_constraints{end+1, 1} = @(k) cellfun(@(C) slice(C, k), grad_dq_ub);
                    end
                    
                    dq_lb_int = interval(dq_lb{i, 1}{j, 1});
                    if ~(dq_lb_int.sup < 0)
                        fprintf('ADDED LOWER BOUND JOINT VELOCITY CONSTRAINT ON JOINT %d AT TIME %d \n', j, i);
                        P.constraints{end+1, 1} = @(k) slice(dq_lb{i, 1}{j, 1}, k);
                        grad_dq_lb = grad(dq_lb{i, 1}{j, 1}, P.jrs_info.n_q);
                        P.grad_constraints{end+1, 1} = @(k) cellfun(@(C) slice(C, k), grad_dq_lb);
                    end
                end
            end

            % more constraints to add:
            % 1) orientation constraints
            % 2) self intersection constraints
            % 3) joint workspace position constraints?
        end

        function [h_obs, grad_h_obs] = evaluate_obs_constraint(P, c, grad_c, k) 
            % made a separate function to handle the obstacle constraints,
            % because the gradient requires knowing the index causing
            % the max of the constraints
            [h_obs, max_idx] = max(c(k));
            h_obs = -h_obs;
            grad_eval = grad_c(k);
            grad_h_obs = zeros(length(k), 1);
            for i = 1:length(k)
                grad_h_obs(i, 1) = -grad_eval{i}(max_idx, :);
            end
        end

        function [h_obs, grad_h_obs] = evaluate_smooth_obs_constraint(P, c, grad_c, k, lambda) 
            % for evaluating smooth obs constraints, where lambda is introduced
            % as extra decision variables to avoid taking a max.

            % evaluate constraint: (A*FO - b)'*lambda
            h_obs = -c(k)'*lambda;

            % evaluate gradient w.r.t. k... 
            % grad_c(k) gives n_k x 1 cell, each containing
            % an N x 1 vector, where N is the number of rows of A.
            % take the dot product of each cell with lambda
            grad_eval = grad_c(k);
            grad_h_obs = zeros(length(k), 1);
            for i = 1:length(k)
                grad_h_obs(i, 1) = -grad_eval{i}'*lambda;
            end

            % evaluate gradient w.r.t. lambda...
            % this is just (A*FO-b)
            grad_h_obs = [grad_h_obs; -c(k)];
        end

        
        function [k_opt, trajopt_failed] = trajopt(P, q_0, q_dot_0, q_ddot_0, q_des)
            % use fmincon to optimize the cost subject to constraints
            P.vdisp('Running trajopt', 3)
            
            P.trajopt_start_tic = tic ;
            n_k = P.jrs_info.n_k;
            if P.smooth_obstacle_constraints_flag && ~isempty(P.smooth_obs_lambda_index)
                cost_func = @(x) P.eval_cost(x(1:n_k), q_0, q_dot_0, q_ddot_0, q_des);
                constraint_func = @(x) P.eval_smooth_constraint(x(1:n_k), x(n_k+1:end));
                lb_k = -ones(n_k, 1);
                ub_k = ones(n_k, 1);
                lb_lambda = zeros(P.smooth_obs_lambda_index{end}(end), 1);
                ub_lambda = ones(P.smooth_obs_lambda_index{end}(end), 1);
                lb = [lb_k; lb_lambda];
                ub = [ub_k; ub_lambda];
            else
                cost_func = @(k) P.eval_cost(k, q_0, q_dot_0, q_ddot_0, q_des);
                constraint_func = @(k) P.eval_constraint(k);
                lb = -ones(n_k, 1);
                ub = ones(n_k, 1);
            end

            initial_guess = rand_range(lb, ub);
           
            options = optimoptions('fmincon','SpecifyConstraintGradient',true);
%             options = optimoptions('fmincon','SpecifyConstraintGradient',true, 'CheckGradients', true);
            [k_opt, ~, exitflag, ~] = fmincon(cost_func, initial_guess, [], [], [], [], lb, ub, constraint_func, options) ;
            
            trajopt_failed = exitflag <= 0 ;
        end
        
        function [cost] = eval_cost(P, k, q_0, q_dot_0, q_ddot_0, q_des)
            if P.use_q_plan_for_cost
                q_plan = P.desired_trajectory(q_0, q_dot_0, q_ddot_0, P.t_plan, k);
                cost = sum((q_plan - q_des).^2);
            else
                % use final position as cost
                q_stop = P.desired_trajectory(q_0, q_dot_0, q_ddot_0, P.t_stop, k);
                cost = sum((q_stop - q_des).^2);
            end
        end

        function [h, heq, grad_h, grad_heq] = eval_constraint(P, k)
            n_obs_c = length(P.obs_constraints);
            n_c = length(P.constraints);

            h = zeros(n_c + n_obs_c, 1);
            grad_h = zeros(length(k), n_c + n_obs_c);
            
            for i = 1:n_obs_c
                [h_i, grad_h_i] = P.obs_constraints{i}(k);
                h(i) = h_i;
                grad_h(:, i) = grad_h_i;
            end

            for i = 1:n_c
                h(i + n_obs_c) = P.constraints{i}(k);
                grad_h(:, i + n_obs_c) = P.grad_constraints{i}(k);
            end

            grad_heq = [];
            heq = [];
        end

        function [h, heq, grad_h, grad_heq] = eval_smooth_constraint(P, k, lambda)

            n_obs_c = length(P.smooth_obs_constraints);
            n_k = length(k);
            n_lambda = length(lambda);

            % max_lambda_index = P.smooth_obs_lambda_index{end}(end);

            % number of constraints: 
            % - 1 obstacle avoidance per obstacle
            % - 1 sum lambdas for obstacle = 1
            % - n_lambda lambda \in {0, 1}

            h = zeros(2*n_obs_c, 1);
            grad_h = zeros(n_k + n_lambda, 2*n_obs_c);

            for i = 1:n_obs_c
                lambda_idx = P.smooth_obs_lambda_index{i};
                lambda_i = lambda(lambda_idx, 1);% pull out correct lambdas!!
                [h_i, grad_h_i] = P.smooth_obs_constraints{i}(k, lambda_i);

                % obs avoidance constraints:
                h(i, 1) = h_i;
                grad_h(1:n_k, i) = grad_h_i(1:n_k, 1);
                grad_h(n_k + lambda_idx, i) = grad_h_i(n_k+1:end, 1);

                % sum lambdas for this obstacle constraint >= 1
                % sum lambdas for this obstacle constraint <= 1?
%                 h(n_obs_c + i, 1) = 1 - sum(lambda_i, 1); 
%                 grad_h(n_k + lambda_idx, n_obs_c + i) = -ones(length(lambda_i), 1);

                % from Borrelli paper
                h(n_obs_c + i, 1) = norm(P.smooth_obs_constraints_A{i}'*lambda_i, 2) - 1;
                % implement gradient here!!!
                A_bar = P.smooth_obs_constraints_A{i}*P.smooth_obs_constraints_A{i}';
                grad_h(n_k + lambda_idx, n_obs_c + i) = 0.5*(lambda_i'*A_bar*lambda_i)^(-0.5)*2*A_bar*lambda_i;
                
            end

            % lambda \in {0, 1} for each lambda
%             heq = lambda.*(lambda - 1);
%             grad_heq = zeros(n_k + n_lambda, n_lambda);
%             grad_heq(n_k+1:end, :) = diag(2*lambda - 1);

             heq = [];
             grad_heq = [];
        end

        function [q_des, qd_des, qdd_des] = desired_trajectory(P, q_0, q_dot_0, q_ddot_0, t, k)
            % at a given time t and traj. param k value, return
            % the desired position, velocity, and acceleration.
            switch P.traj_type
            case 'orig'
                t_plan = P.t_plan;
                t_stop = P.t_stop;
                k_scaled = P.jrs_info.c_k_orig + P.jrs_info.g_k_orig.*k;
                
                if ~isnan(k)
                    if t <= t_plan
                        % compute first half of trajectory
                        q_des = q_0 + q_dot_0*t + (1/2)*k_scaled*t^2;
                        qd_des = q_dot_0 + k_scaled*t;
                        qdd_des = k_scaled;
                    else
                        % compute trajectory at peak
                        q_peak = q_0 + q_dot_0*t_plan + (1/2)*k_scaled*t_plan^2;
                        qd_peak = q_dot_0 + k_scaled*t_plan;

                        % compute second half of trajectory
                        q_des = q_peak + qd_peak*(t-t_plan) + (1/2)*((0 - qd_peak)/(t_stop - t_plan))*(t-t_plan)^2;
                        qd_des = qd_peak + ((0 - qd_peak)/(t_stop - t_plan))*(t-t_plan);
                        qdd_des = (0 - qd_peak)/(t_stop - t_plan);
                    end
                else
                    % bring the trajectory to a stop in t_plan seconds
                    % trajectory peaks at q_0
                    q_peak = q_0;
                    qd_peak = q_dot_0;
                    
                    if t <= t_plan && ~all(q_dot_0 == 0) % we're braking!
                        q_des = q_peak + qd_peak*t + (1/2)*((0 - qd_peak)/t_plan)*t^2;
                        qd_des = qd_peak + ((0 - qd_peak)/t_plan)*t;
                        qdd_des = (0 - qd_peak)/t_plan;
                    else % we should already be stopped, maintain that.
                        q_des = q_peak;
                        qd_des = zeros(size(q_dot_0));
                        qdd_des = zeros(size(q_0));
                    end
                end
            case 'bernstein'
                % assuming K = [-1, 1] corresponds to final position for now!!
                n_q = length(q_0);
                if ~isnan(k)
                    q1 = q_0 + P.jrs_info.c_k_bernstein + P.jrs_info.g_k_bernstein.*k;
                    for j = 1:n_q
                        beta{j} = match_deg5_bernstein_coefficients({q_0(j); q_dot_0(j); q_ddot_0(j); q1(j); 0; 0});
                        alpha{j} = bernstein_to_poly(beta{j}, 5);
                    end
                    q_des = zeros(length(q_0), 1);
                    qd_des = zeros(length(q_0), 1);
                    qdd_des = zeros(length(q_0), 1);
                    for j = 1:n_q
                        for coeff_idx = 0:5
                            q_des(j) = q_des(j) + alpha{j}{coeff_idx+1}*t^coeff_idx;
                            if coeff_idx > 0
                                qd_des(j) = qd_des(j) + coeff_idx*alpha{j}{coeff_idx+1}*t^(coeff_idx-1);
                            end
                            if coeff_idx > 1
                                qdd_des(j) = qdd_des(j) + (coeff_idx)*(coeff_idx-1)*alpha{j}{coeff_idx+1}*t^(coeff_idx-2);
                            end
                        end
                    end
                else
                    % bring the trajectory to a stop using previous trajectory...
                    t_plan = P.t_plan;
                    if t <= t_plan && norm(q_dot_0) > 1e-8 && norm(q_dot_0) > 1e-8
                        % just plug into previous trajectory, but shift time forward by t_plan.
                        [q_des, qd_des, qdd_des] = P.info.desired_trajectory{end - 1}(t + t_plan);
                    else % we should already be stopped, maintain that.
                        q_des = q_0;
                        qd_des = zeros(n_q, 1);
                        qdd_des = zeros(n_q, 1);
                    end
                end
            otherwise
                error('trajectory type not recognized');
            end
        end
        
    end
end

