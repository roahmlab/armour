function [u, f, n] = rnea(q, qd, q_aux_d, qdd, use_gravity, robot_params)
    %% Robot parameters
    % link masses
    mass = robot_params.mass;
    
    % center of mass for each link
    com = robot_params.com;
    
    % inertia wrt to center of mass frame
    I = robot_params.I;
    
    % number of joints
    num_joints = robot_params.num_joints;
    
    % number of active joints
    num_q = robot_params.num_q;
    
    % use interval arithmetic
    use_interval = strcmp(robot_params.set_type, 'interval');
    
    % fixed transforms
    T0 = robot_params.T0;
    joint_axes = robot_params.joint_axes;
    joint_types = robot_params.joint_types;
    
    %% set inputs
    joint_pos = q(:);
    joint_vel = qd(:);
    joint_acc = qdd(:);
    joint_vel_aux = q_aux_d(:);

    %% setup reference frames
    % rotation axis of base frame
    z0 = [0;0;1];

    % orientation of frame i with respect to frame i-1
    R = repmat(eye(3), [1, 1, num_joints+1]);

    % position of the origin of frame i with respect to frame i-1
    P = zeros(3, num_joints+1);

    % orientation of frame i-1 with respect to frame i
    R_t = repmat(eye(3), [1, 1, num_joints]);
    
    % Frame {i} axis of rotation expressed in Frame {i}
    z = zeros(3, num_joints);

    % calculate frame-to-frame transformations based on DH table
    for i = 1:num_joints
        switch joint_types{i}
            case 'revolute'
                % orientation of Frame {i} with respect to Frame {i-1}
                dim = find(joint_axes(:,i) ~=0);
                if dim == 1
                    R(:,:,i) = T0(1:3, 1:3, i) * rx(joint_pos(robot_params.q_index(i)));
                elseif dim == 2
                    R(:,:,i) = T0(1:3, 1:3, i) * ry(joint_pos(robot_params.q_index(i)));
                else
                    R(:,:,i) = T0(1:3, 1:3, i) * rz(joint_pos(robot_params.q_index(i)));
                end
            case 'fixed'
                R(:, :, i) = T0(1:3, 1:3, i);
            otherwise
                error('Joint type not implemented');
        end
        
        R_t(:,:,i) = R(:,:,i)'; % line 7
 
        % position of Frame {i} with respect to Frame {i-1}
        P(:,i) = T0(1:3, 4, i);
        
        % orientation of joint i axis of rotation with respect to Frame {i}
        z(:,i) = joint_axes(:,i);
    end
 
    % get transform to end-effector
    if robot_params.num_bodies > robot_params.num_joints
        R(:, :, end) = T0(1:3, 1:3, end);
        P(:, end) = T0(1:3, 4, end);
    end
    
    %% INITIALIZE
    % base link/frame
    w0 = zeros(3,1);
    w0dot = zeros(3,1);
    linear_acc0 = zeros(3,1);
    w0_aux = zeros(3,1); % auxilliary

    % set gravity
    if use_gravity == true
        linear_acc0(:,1) = -robot_params.gravity;
    end

    % angular velocity/acceleration
    w = zeros(3, num_joints);
    wdot = zeros(3, num_joints);
    w_aux = zeros(3, num_joints);

    % linear acceleration of frame
    linear_acc = zeros(3, num_joints);
    
    if ~use_interval
        
        % linear acceleration of com
        linear_acc_com = zeros(3, num_joints);

        % link forces/torques
        F = zeros(3, num_joints);
        N = zeros(3, num_joints);
        
        % intialize f, n, u
        f = zeros(3, num_joints + 1);
        n = zeros(3, num_joints + 1);
        u = zeros(num_q,1);
    else
        % linear acceleration of com
        linear_acc_com = interval(zeros(3, num_joints), zeros(3, num_joints));

        % link forces/torques
        F = interval(zeros(3, num_joints), zeros(3, num_joints));
        N = interval(zeros(3, num_joints), zeros(3, num_joints));
        
        % intialize f, n, u
        f = interval(zeros(3, num_joints + 1), zeros(3, num_joints + 1));
        n = interval(zeros(3, num_joints + 1), zeros(3, num_joints + 1));
        u = interval(zeros(num_q, 1), zeros(num_q,1));
    end
    
    %% RNEA forward recursion
    for i=1:num_joints
        % MATLAB has no zero indexing, soo....
        if i ==1
            switch joint_types{i}
                case 'revolute'
                    % (6.45) angular velocity
                    w(:,i) = R_t(:,:,i) * w0 + joint_vel(robot_params.q_index(i))*z(:,i); % line 13

                    % auxillary angular velocity
                    w_aux(:,i) = R_t(:,:,i) * w0_aux + joint_vel_aux(robot_params.q_index(i))*z(:,i); % line 13

                    % (6.46) angular acceleration
                    wdot(:,i) = R_t(:,:,i) * w0dot ... % line 15
                                + cross(R_t(:,:,i) * w0_aux, joint_vel(robot_params.q_index(i))*z(:,i))...
                                + joint_acc(robot_params.q_index(i))*z(:,i);

                    % (6.47) linear acceleration        
                    linear_acc(:,i) = R_t(:,:,i)*(linear_acc0 ...
                                        + cross(w0dot, P(:,i)) ... % line 16 (TYPO IN PAPER)
                                        + cross(w0, cross(w0, P(:,i))));      
                case 'fixed'
                    % (6.45) angular velocity
                    w(:,i) = R_t(:,:,i) * w0; % line 13

                    % auxillary angular velocity
                    w_aux(:,i) = R_t(:,:,i) * w0_aux; % line 13

                    % (6.46) angular acceleration
                    wdot(:,i) = R_t(:,:,i) * w0dot; % line 15

                    % (6.47) linear acceleration        
                    linear_acc(:,i) = R_t(:,:,i)*(linear_acc0 ...
                                        + cross(w0dot, P(:,i)) ... % line 16 (TYPO IN PAPER)
                                        + cross(w0, cross(w0, P(:,i))));  
                otherwise
                    error('Joint type not implemented');
            end
                            
        else
            switch joint_types{i}
                case 'revolute'
                    % (6.45) angular velocity
                    w(:,i) = R_t(:,:,i) * w(:,i-1) + joint_vel(robot_params.q_index(i))*z(:,i); % line 13

                    % auxillary angular velocity
                    w_aux(:,i) = R_t(:,:,i) * w_aux(:,i-1) + joint_vel_aux(robot_params.q_index(i))*z(:,i); % line 14

                    % (6.46) angular acceleration
                    wdot(:,i) = R_t(:,:,i) * wdot(:,i-1) ... % line 15
                                        + cross(R_t(:,:,i) * w_aux(:,i-1), joint_vel(robot_params.q_index(i))*z(:,i))...
                                        + joint_acc(robot_params.q_index(i))*z(:,i);

                    % (6.47) linear acceleration
                    linear_acc(:,i) = R_t(:,:,i)*(linear_acc(:,i-1) ...
                                            + cross(wdot(:,i-1), P(:,i)) ... % line 16 (TYPO IN PAPER)
                                            + cross(w(:,i-1), cross(w_aux(:,i-1),P(:,i))));
                case 'fixed'
                    % (6.45) angular velocity
                    w(:,i) = R_t(:,:,i) * w(:,i-1); % line 13

                    % auxillar angular velocity
                    w_aux(:,i) = R_t(:,:,i) * w_aux(:,i-1); % line 14

                    % (6.46) angular acceleration
                    wdot(:,i) = R_t(:,:,i) * wdot(:,i-1); % line 15

                    % (6.47) linear acceleration
                    linear_acc(:,i) = R_t(:,:,i)*(linear_acc(:,i-1) ...
                                            + cross(wdot(:,i-1), P(:,i)) ... % line 16 (TYPO IN PAPER)
                                            + cross(w(:,i-1), cross(w_aux(:,i-1),P(:,i))));
                otherwise
                    error('Joint type not implemented');
            end

        end  

        % (6.48) linear acceleration of CoM auxilliary
        linear_acc_com(:,i) = linear_acc(:,i) ... % line 23 (modified for standard RNEA)
                        + cross(wdot(:,i),com(:,i)) .... 
                        + cross(w(:,i),cross(w_aux(:,i),com(:,i)));                

        % (6.49) calculate forces
        F(:,i) = mass(:,i) * linear_acc_com(:,i); % line 27

        % (6.50) calculate torques
        N(:,i) = I{i}*wdot(:,i) ... % calculated in line 29
                 + cross(w_aux(:,i), (I{i}*w(:,i)));

    end

    %% RNEA reverse recursion
    for i = num_joints:-1:1
        % (6.51)
        f(:,i) = R(:,:,i+1) * f(:,i+1) + F(:,i); % line 28

        % (6.52)
        n(:,i) = N(:,i) ...
               + R(:,:,i+1) * n(:,i+1) ... % line 29
               + cross(com(:,i), F(:,i)) ... % P(:,i) might not be right
               + cross(P(:,i+1), R(:,:,i+1)*f(:,i+1)); % line 29 (TYPO IN PAPER)
    end

    % calculate joint torques
    for i = 1:num_joints
        switch joint_types{i}
            case 'revolute'
                % (6.53)
                u(robot_params.q_index(i),1) = n(:,i)' * z(:,i); % line 31
            case 'fixed'
                continue;
            otherwise
                error('Joint type not implemented');
        end
    end
end
