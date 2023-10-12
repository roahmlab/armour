close all; clear; clc;

%% user parameters
use_robust_input = true;

dimension = 3 ;
verbosity = 10;

%%% for planner
traj_type = 'orig'; % pick 'orig' or 'bernstein'
allow_replan_errors = true ;
first_iter_pause_flag = false;
use_q_plan_for_cost = false; % otherwise use q_stop (q at final time)
input_constraints_flag = true;

%%% for agent
agent_urdf = 'kinova_without_gripper.urdf';
model_uncertainties = 0.0 : 0.05 : 0.3;
add_uncertainty_to = 'all'; % choose 'all', 'link', or 'none'
links_with_uncertainty = {}; % if add_uncertainty_to = 'link', specify links here.


agent_move_mode = 'integrator' ; % pick 'direct' or 'integrator'
use_CAD_flag = false;
add_measurement_noise_ = false;
measurement_noise_size_ = 0;

%%% for LLC
LLC_V_max = 1e-2;
use_true_params_for_robust = false;
if_use_mex_controller = true;

for model_uncertainty = model_uncertainties
%% robot params:
    uncertain_mass_range = [1 - model_uncertainty, 1 + model_uncertainty];
    
    robot = importrobot(agent_urdf);
    robot.DataFormat = 'col';
    robot.Gravity = [0 0 -9.81];
    params = load_robot_params(robot, ...
                               'add_uncertainty_to', add_uncertainty_to, ...
                               'links_with_uncertainty', links_with_uncertainty,...
                               'uncertain_mass_range', uncertain_mass_range);
    joint_speed_limits = [-1.3963, -1.3963, -1.3963, -1.3963, -1.2218, -1.2218, -1.2218;
                           1.3963,  1.3963,  1.3963,  1.3963,  1.2218,  1.2218,  1.2218]; % matlab doesn't import these from urdf so hard code into class
    joint_input_limits = [-56.7, -56.7, -56.7, -56.7, -29.4, -29.4, -29.4;
                           56.7,  56.7,  56.7,  56.7,  29.4,  29.4,  29.4]; % matlab doesn't import these from urdf so hard code into class
    transmision_inertia = [8.02999999999999936 11.99620246153036440 9.00254278617515169 11.58064393167063599 8.46650409179141228 8.85370693737424297 8.85873036646853151]; % matlab doesn't import these from urdf so hard code into class
    M_min_eigenvalue = 5.095620491878957; % matlab doesn't import these from urdf so hard code into class
    eps = sqrt(2 * LLC_V_max / M_min_eigenvalue);
    
    % create arm agent
    A = uarmtd_agent(robot, params,...
                     'verbose', verbosity,...
                     'animation_set_axes_flag', 0,... 
                     'animation_set_view_flag', 0,...
                     'move_mode', agent_move_mode,...
                     'use_CAD_flag', use_CAD_flag,...
                     'joint_speed_limits', joint_speed_limits, ...
                     'joint_input_limits', joint_input_limits, ...
                     'add_measurement_noise_', add_measurement_noise_, ...
                     'measurement_noise_size_', measurement_noise_size_,...
                     'M_min_eigenvalue', M_min_eigenvalue, ...
                     'transmision_inertia', transmision_inertia);
    
    %% automated from here
    robust_out_armour_final = [];
    robust_out_althoff_final = [];
    
    for sam = 1:100
        fprintf('%.2f %d\n', model_uncertainty, sam);
    
        %% initial conditions
        q0 = pi * rand(A.n_inputs, 1) - pi / 2;
        qd0 = pi * rand(A.n_inputs, 1) - pi / 2;
        
        z0 = zeros(A.n_states,1);
        z0(A.joint_state_indices) = q0;
        z0(A.joint_speed_indices) = qd0;
        
        %% trajectory params
        q_diff = rand(A.n_inputs, 1) * 2 - 1; q_diff = 0.025 * pi * q_diff / norm(q_diff);
        qd_diff = rand(A.n_inputs, 1) * 2 - 1; qd_diff = 0.05 * pi * qd_diff / norm(qd_diff);
        traj_q0 = q0 + q_diff;
        traj_qd0 = qd0 + qd_diff;
        traj_qdd0 = zeros(A.n_inputs, 1);
        T = 2.5;
        
        planner_info.desired_trajectory = {@(t)desired_trajectory(traj_q0, traj_qd0, traj_qdd0, t, T)};
        
        %% run simulation amd process data (armour)
        A.LLC = uarmtd_robust_CBF_LLC('verbose', verbosity, ...
                                      'use_true_params_for_robust', use_true_params_for_robust, ...
                                      'V_max', LLC_V_max, ...
                                      'if_use_mex_controller', if_use_mex_controller);
        
        [tout_armour, zout_armour] = A.integrator(@(t,z) A.dynamics(t,z,planner_info), [0, T], z0);
        
        q_armour = zout_armour(A.joint_state_indices,:);
        qd_armour = zout_armour(A.joint_speed_indices,:);
        
        uout_armour = zeros(A.n_inputs, size(tout_armour, 2));
        nominal_out_armour = zeros(size(uout_armour));
        robust_out_armour = zeros(size(uout_armour));
        
        % store approximate inputs at each time:
        for j = 1:length(tout_armour)
            t = tout_armour(j);
            z_meas = zout_armour(:,j);
            [uout_armour(:, j), nominal_out_armour(:, j), robust_out_armour(:, j)] = A.LLC.get_control_inputs(A, t, z_meas, planner_info);
        end
        
        robust_out_armour_final = [robust_out_armour_final, robust_out_armour];
    
        %% run simulation amd process data (althoff)
        A.LLC = uarmtd_robust_CBF_LLC('verbose', verbosity, ...
                                      'use_true_params_for_robust', use_true_params_for_robust, ...
                                      'V_max', LLC_V_max, ...
                                      'if_use_mex_controller', if_use_mex_controller, ...
                                      'if_use_armour_robust_controller', false);
        
        [tout_althoff, zout_althoff] = A.integrator(@(t,z) A.dynamics(t,z,planner_info), [0, T], z0);
        
        q_althoff = zout_althoff(A.joint_state_indices,:);
        qd_althoff = zout_althoff(A.joint_speed_indices,:);
        
        uout_althoff = zeros(A.n_inputs, size(tout_armour, 2));
        nominal_out_althoff = zeros(size(uout_armour));
        robust_out_althoff = zeros(size(uout_armour));
        
        % store approximate inputs at each time:
        for j = 1:length(tout_althoff)
            t = tout_althoff(j);
            z_meas = zout_althoff(:,j);
            [uout_althoff(:, j), nominal_out_althoff(:, j), robust_out_althoff(:, j)] = A.LLC.get_control_inputs(A, t, z_meas, planner_info);
        end
    
        robust_out_althoff_final = [robust_out_althoff_final, robust_out_althoff];
    end
    
    %% comparison
    % [q_des, qd_des] = desired_trajectory(traj_q0, traj_qd0, traj_qdd0, tout_armour, T);
    % 
    % figure;
    % for i = 1:7
    %     subplot(3,3,i); hold on;
    %     plot(tout_althoff, robust_out_althoff(i,:), 'r');
    %     plot(tout_armour, robust_out_armour(i,:), 'b');
    %     xlabel('time (sec)');
    %     ylabel('robust input (N*m)');
    %     if i == 1
    %         legend('Althoff','Armour');
    %     end
    %     title(['joint', num2str(i)]);
    % end
    % sgtitle(['robust input (model uncertainty ', num2str(100 * model_uncertainty), ' %)']);
    % 
    % figure;
    % for i = 1:7
    %     subplot(3,3,i); hold on;
    %     plot(tout_althoff, rad2deg(wrapToPi(q_althoff(i,:) - q_des(i,:))), 'r');
    %     plot(tout_armour, rad2deg(wrapToPi(q_armour(i,:) - q_des(i,:))), 'b');
    %     plot(tout_armour, -rad2deg(eps / A.LLC.Kr) * ones(1,length(tout_armour)), 'k');
    %     plot(tout_armour, rad2deg(eps / A.LLC.Kr) * ones(1,length(tout_armour)), 'k');
    %     xlabel('time (sec)');
    %     ylabel('tracking error (degree)');
    %     if i == 1
    %         legend('Althoff','Armour','ultimate bound','ultimate bound');
    %     end
    %     title(['joint ', num2str(i)]);
    % end
    % sgtitle(['tracking error (model uncertainty ', num2str(100 * model_uncertainty), ' %)']);
    
    save(['ultimate_robust_input_', num2str(100 * model_uncertainty), '.mat'], 'robust_out_armour_final', 'robust_out_althoff_final');
end

%% helper functions
function [q_des, qd_des, qdd_des] = desired_trajectory(q0, qd0, qdd0, t, T)
    beta = match_deg5_bernstein_coefficients({q0, qd0, qdd0, 0, 0, 0}, T);

    [B, dB, ddB] = Bezier_kernel_deg5(t / T);

    q_des = zeros(length(q0), length(t));
    qd_des = zeros(length(qd0), length(t));
    qdd_des = zeros(length(qdd0), length(t));

    for i = 1:length(beta)
        q_des = q_des + beta{i} * B(:,i)';
        qd_des = qd_des + beta{i} * dB(:,i)';
        qdd_des = qdd_des + beta{i} * ddB(:,i)';
    end

    qd_des = qd_des / T;
    qdd_des = qdd_des / T / T;
end



















