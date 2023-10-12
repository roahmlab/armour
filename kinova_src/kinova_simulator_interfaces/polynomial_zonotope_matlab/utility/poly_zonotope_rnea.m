function [u, f, n] = poly_zonotope_rnea(R_in, R_t_in, qd, q_aux_d, qdd, use_gravity, robot_params)
    % outputs: u - nx1 cell of tortatotopes for n joints

	%% Robot parameters
	% max number of generators to keep around...
	% order * dimension = max # of generators
	zono_order = robot_params.zono_order; % <-- change this to tradeoff speed vs conservatism
    
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
    
    % fixed transforms
    T0 = robot_params.T0;
    joint_axes = robot_params.joint_axes;
    joint_types = robot_params.joint_types;
    
    %% set inputs
%     joint_pos = q;
    joint_vel = qd;
    joint_acc = qdd;
    joint_vel_aux = q_aux_d;
    
    %% reference frames
    % rotation axis of base frame
    z0 = [0;0;1];

    % position of the origin of frame i with respect to frame i-1
    P = robot_params.P;

    % Frame {i} axis of rotation expressed in Frame {i}
    z = robot_params.joint_axes;
    
    % orientation of frame i with respect to frame i-1 (and transpose)
    % generate matPolyZonotopes for R from q
%     [R, R_t] = get_transforms_from_jrs(q, z);
    for i = 1:num_joints
        switch joint_types{i}
            case 'revolute'
                % [R{i, 1}, R_t{i, 1}] = get_pz_rotations_from_q(joint_pos{robot_params.q_index(i)}, z(:, i));
                R{i, 1} = matPolyZonotope_ROAHM(T0(1:3, 1:3, i)) * R_in{robot_params.q_index(i), 1};
                R_t{i, 1} = R_t_in{robot_params.q_index(i), 1} * matPolyZonotope_ROAHM(T0(1:3, 1:3, i)');
            case 'fixed'
                R{i, 1} = matPolyZonotope_ROAHM(T0(1:3, 1:3, i));
                R_t{i, 1} = matPolyZonotope_ROAHM(R{i, 1}.C');
            otherwise
                error('Joint type not implemented');
        end
    end
    
    % get transform to end-effector
%     if robot_params.num_bodies > robot_params.num_joints
%         R{end+1, 1} = T0(1:3, 1:3, end);
%         R_t{end+1, 1} = R{end, 1}';
%         P(:, end+1) = T0(1:3, 4, end);
%     end

    %% INITIALIZE
    % base link/frame 
    w0 = zeros(3,1);
    w0dot = zeros(3,1);
    linear_acc0 = zeros(3,1);
    w0_aux = zeros(3,1); % auxilliary

    % set gravity
    if use_gravity == true
        linear_acc0 = -robot_params.gravity';
    end

    %% RNEA forward recursion
    for i=1:num_joints
        % MATLAB has no zero indexing, soo....
        if i ==1
            switch joint_types{i}
                case 'revolute'
                    % (6.45) angular velocity
                    w{i, 1} = R_t{i} * w0 + joint_vel{robot_params.q_index(i)}*z(:, i); % line 13

                    % auxillary angular velocity
                    w_aux{i, 1} = R_t{i} * w0_aux + joint_vel_aux{robot_params.q_index(i)}*z(:, i); % line 13

                    % (6.46) angular acceleration
                    wdot{i, 1} = R_t{i}* w0dot ... % line 15
                                + cross(R_t{i} * w0_aux, joint_vel{robot_params.q_index(i)}*z(:, i))...
                                + joint_acc{robot_params.q_index(i)}*z(:, i);

                    % (6.47) linear acceleration        
                    linear_acc{i, 1} = R_t{i}*(linear_acc0 ...
                                        + cross(w0dot, P(:, i)) ... % line 16 (TYPO IN PAPER)
                                        + cross(w0, cross(w0, P(:, i))));

                    % reduction step to manage number of generators:
                    w{i, 1} = reduce(w{i, 1}, 'girard', zono_order);
                    w_aux{i, 1} = reduce(w_aux{i, 1}, 'girard', zono_order);
                    wdot{i, 1} = reduce(wdot{i, 1}, 'girard', zono_order);
                    linear_acc{i, 1} = reduce(linear_acc{i, 1}, 'girard', zono_order);  
                case 'fixed'
                    % (6.45) angular velocity
                    w{i, 1} = R_t{i} * w0; % line 13

                    % auxillary angular velocity
                    w_aux{i, 1} = R_t{i} * w0_aux; % line 13

                    % (6.46) angular acceleration
                    wdot{i, 1} = R_t{i}* w0dot; % line 15

                    % (6.47) linear acceleration        
                    linear_acc{i, 1} = R_t{i}*(linear_acc0 ...
                                        + cross(w0dot, P(:, i)) ... % line 16 (TYPO IN PAPER)
                                        + cross(w0, cross(w0, P(:, i))));
                                    
                    % reduction step to manage number of generators:
                    w{i, 1} = reduce(w{i, 1}, 'girard', zono_order);
                    w_aux{i, 1} = reduce(w_aux{i, 1}, 'girard', zono_order);
                    wdot{i, 1} = reduce(wdot{i, 1}, 'girard', zono_order);
                    linear_acc{i, 1} = reduce(linear_acc{i, 1}, 'girard', zono_order);  
                otherwise
                    error('Joint type not implemented');
            end
        else
            switch joint_types{i}
                case 'revolute'
                    % (6.45) angular velocity
                    w{i, 1} = R_t{i} * w{i-1, 1} + joint_vel{robot_params.q_index(i)}*z(:, i); % line 13
                    w{i, 1} = reduce(w{i, 1}, 'girard', zono_order);

                    % auxillary angular velocity
                    w_aux{i, 1} = R_t{i} * w_aux{i-1 ,1} + joint_vel_aux{robot_params.q_index(i)}*z(:, i); % line 14
                    w_aux{i, 1} = reduce(w_aux{i, 1}, 'girard', zono_order);

                    % (6.46) angular acceleration
                    prod1 = R_t{i} * w_aux{i-1, 1}; 
                    prod1 = reduce(prod1, 'girard', zono_order); % <-- these intermediate products should help speed up computation
                    prod2 = joint_vel{robot_params.q_index(i)}*z(:, i);
                    prod2 = reduce(prod2, 'girard', zono_order);
                    wdot{i, 1} = R_t{i, 1} * wdot{i-1, 1}... % line 15
                                        + cross(prod1, prod2)...
                                        + joint_acc{robot_params.q_index(i)}*z(:, i);
                    wdot{i, 1} = reduce(wdot{i, 1}, 'girard', zono_order);

                    % (6.47) linear acceleration
                    linear_acc{i, 1} = R_t{i, 1}*(linear_acc{i-1, 1} ...
                                            + cross(wdot{i-1, 1}, P(:, i)) ... % line 16 (TYPO IN PAPER)
                                            + cross(w{i-1, 1}, cross(w_aux{i-1, 1},P(:, i))));
                    linear_acc{i, 1} = reduce(linear_acc{i, 1}, 'girard', zono_order);
                case 'fixed'
                    % (6.45) angular velocity
                    w{i, 1} = R_t{i} * w{i-1, 1}; % line 13
                    w{i, 1} = reduce(w{i, 1}, 'girard', zono_order);

                    % auxillary angular velocity
                    w_aux{i, 1} = R_t{i} * w_aux{i-1 ,1}; % line 14
                    w_aux{i, 1} = reduce(w_aux{i, 1}, 'girard', zono_order);

                    % (6.46) angular acceleration
                    wdot{i, 1} = R_t{i, 1} * wdot{i-1, 1}; % line 15
                    wdot{i, 1} = reduce(wdot{i, 1}, 'girard', zono_order);

                    % (6.47) linear acceleration
                    linear_acc{i, 1} = R_t{i, 1}*(linear_acc{i-1, 1} ...
                                            + cross(wdot{i-1, 1}, P(:, i)) ... % line 16 (TYPO IN PAPER)
                                            + cross(w{i-1, 1}, cross(w_aux{i-1, 1},P(:, i))));
                    linear_acc{i, 1} = reduce(linear_acc{i, 1}, 'girard', zono_order);
                otherwise
                    error('Joint type not implemented');
            end
            
        end  

        % (6.48) linear acceleration of CoM auxilliary
        prod0 = cross(wdot{i, 1},com{i, 1});
        prod0 = reduce(prod0, 'girard', zono_order);
        prod1 = cross(w_aux{i, 1},com{i, 1});
        prod1 = reduce(prod1, 'girard', zono_order);
        linear_acc_com{i, 1} = linear_acc{i, 1} ... % line 23 (modified for standard RNEA)
                        + prod0 .... 
                        + cross(w{i, 1}, prod1);   
        linear_acc_com{i, 1} = reduce(linear_acc_com{i, 1}, 'girard', zono_order);

        % (6.49) calculate forces
        F{i, 1} = mass{i, 1} * linear_acc_com{i, 1}; % line 27
        F{i, 1} = reduce(F{i, 1}, 'girard', zono_order);

        % (6.50) calculate torques
        prod0 = I{i}*wdot{i, 1};
        prod0 = reduce(prod0, 'girard', zono_order);
        prod1 = (I{i}*w{i, 1});
        prod1 = reduce(prod1, 'girard', zono_order);
        N{i ,1} =  prod0 ... % calculated in line 29
                 + cross(w_aux{i, 1}, prod1);
        N{i, 1} = reduce(N{i, 1}, 'girard', zono_order);
    end

    %% RNEA reverse recursion
    f{num_joints+1, 1} = zeros(3, 1);
    n{num_joints+1, 1} = zeros(3, 1);
    R{end+1, 1} = matPolyZonotope_ROAHM(T0(1:3, 1:3, end));
    R_t{end+1, 1} = matPolyZonotope_ROAHM(R{end, 1}.C');
    P(:, end+1) = T0(1:3, 4, end);
    for i = num_joints:-1:1
        % (6.51)
        f{i, 1} = R{i+1} * f{i+1, 1} + F{i, 1}; % line 28
        f{i, 1} = reduce(f{i, 1}, 'girard', zono_order);

        % (6.52)
        prod0 = R{i+1, 1} * n{i+1, 1};
        prod0 = reduce(prod0, 'girard', zono_order);
        prod1 = R{i+1, 1}*f{i+1, 1};
        prod1 = reduce(prod1, 'girard', zono_order);
        n{i, 1} = N{i, 1} ...
               + prod0 ... % line 29
               + cross(com{i, 1}, F{i, 1}) ... % P(:,i) might not be right
               + cross(P(:, i+1), prod1); % line 29 (TYPO IN PAPER)
        n{i, 1} = reduce(n{i, 1}, 'girard', zono_order);
    end

    % calculate joint torques
    for i = 1:num_joints
        switch joint_types{i}
            case 'revolute'
                % (6.53)
                u{robot_params.q_index(i), 1} = z(:, i)'*n{i, 1};
            case 'fixed'
                continue;
            otherwise
                error('Joint typpe not implemented');
        end
    end
end