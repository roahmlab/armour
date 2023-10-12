% function [u, f, n] = rnea(q, qd, q_aux_d, qdd, use_gravity, robot_params)

clear; clc;

robot_name = 'kinova';
agent_urdf = 'kinova_without_gripper.urdf';
robot = importrobot(agent_urdf);
robot.DataFormat = 'col';
robot.Gravity = [0 0 -9.81];
params = load_robot_params(robot);
robot_params = params.nominal;

q = sym('q', [robot_params.num_q, 1], 'real');
qd = sym('qd', [robot_params.num_q, 1], 'real');
q_aux_d = sym('qda', [robot_params.num_q, 1], 'real');
qdd = sym('qdd', [robot_params.num_q, 1], 'real');

%% Robot parameters
% link masses
% mass = robot_params.mass;
mass = sym('mass', size(robot_params.mass), 'real');

% center of mass for each link
% com = robot_params.com;
com = sym('com', size(robot_params.com), 'real');

% inertia wrt to center of mass frame
% I = robot_params.I;
I = cell(length(robot_params.I),1);
for i = 1:length(robot_params.I)
    I{i} = sym(['I', num2str(i), '_'], [3,3], 'real');
end

% number of joints
num_joints = robot_params.num_joints;

% number of active joints
num_q = robot_params.num_q;

% use interval arithmetic
use_interval = strcmp(robot_params.set_type, 'interval');

% fixed transforms
T0 = robot_params.T0;

threshold = 1e-7;

for i = 1:size(T0, 3)
    R = T0(1:3,1:3,i);
    p = T0(1:3,4,i);
    
    R(abs(R) <= threshold) = 0;
    R(abs(R - 1) <= threshold) = 1;
    R(abs(R + 1) <= threshold) = -1;
    p(abs(p) <= threshold) = 0;
    
    T0(1:3,1:3,i) = R;
    T0(1:3,4,i) = p;
end

joint_axes = robot_params.joint_axes;
joint_types = robot_params.joint_types;

%% set inputs
joint_pos = q(:);
joint_vel = qd(:);
joint_acc = qdd(:);
joint_vel_aux = q_aux_d(:);

%% setup reference frames
% rotation axis of base frame
z0 = sym([0;0;1]);

% orientation of frame i with respect to frame i-1
R = sym(repmat(eye(3), [1, 1, num_joints+1]));

% position of the origin of frame i with respect to frame i-1
P = sym(zeros(3, num_joints+1));

% orientation of frame i-1 with respect to frame i
R_t = sym(repmat(eye(3), [1, 1, num_joints]));

% Frame {i} axis of rotation expressed in Frame {i}
z = sym(zeros(3, num_joints));

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

cosq = sym('cq', [robot_params.num_joints, 1], 'real');
sinq = sym('sq', [robot_params.num_joints, 1], 'real');

R = subs(R, cos(q), cosq);
R = subs(R, sin(q), sinq);
R_t = subs(R_t, cos(q), cosq);
R_t = subs(R_t, sin(q), sinq);

%% INITIALIZE
% base link/frame
w0 = sym(zeros(3,1));
w0dot = sym(zeros(3,1));
linear_acc0 = sym(zeros(3,1));
w0_aux = sym(zeros(3,1)); % auxilliary

% set gravity
linear_acc0(:,1) = -robot_params.gravity;

% angular velocity/acceleration
w = sym('w', [3, num_joints]);
wdot = sym('wdot', [3, num_joints]);
w_aux = sym('w_aux', [3, num_joints]);

% linear acceleration of frame
linear_acc = sym('linear_acc', [3, num_joints]);

w_old = sym('w_old', [3, 1]);
wdot_old = sym('wdot_old', [3, 1]);
w_aux_old = sym('w_aux_old', [3, 1]);
linear_acc_old = sym('linear_acc_old', [3, 1]);

% linear acceleration of com
linear_acc_com = sym(zeros(3, num_joints));

% link forces/torques
F = sym(zeros(3, num_joints));
N = sym(zeros(3, num_joints));

% intialize f, n, u
f = sym(zeros(3, num_joints + 1));
n = sym(zeros(3, num_joints + 1));
u = sym(zeros(num_q,1));

%% RNEA forward recursion
filename = [robot_name,'_rnea.cpp'];
fid = fopen(filename, 'w');

for i = 1:7
    fprintf(fid, 'PZsparse& cq%d = traj->cos_q_des(%d, s_ind);\n', i,i-1 );
end
fprintf(fid, '\n');
for i = 1:7
    fprintf(fid, 'PZsparse& sq%d = traj->sin_q_des(%d, s_ind);\n', i,i-1 );
end
fprintf(fid, '\n');
for i = 1:7
    fprintf(fid, 'PZsparse& qd%d = traj->qd_des(%d, s_ind);\n', i,i-1 );
end
fprintf(fid, '\n');
for i = 1:7
    fprintf(fid, 'PZsparse& qda%d = traj->qda_des(%d, s_ind);\n', i,i-1 );
end
fprintf(fid, '\n');
for i = 1:7
    fprintf(fid, 'PZsparse& qdd%d = traj->qdda_des(%d, s_ind);\n', i,i-1 );
end
fprintf(fid, '\n');

for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse w', num2str(j),'(0);']);
end
for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse w_aux', num2str(j),'(0);']);
end
for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse wdot', num2str(j),'(0);']);
end
for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse linear_acc', num2str(j),'(0);']);
end
fprintf(fid, '\n');

for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse w_new', num2str(j),'(0);']);
end
for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse w_aux_new', num2str(j),'(0);']);
end
for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse wdot_new', num2str(j),'(0);']);
end
for j = 1:3
    fprintf(fid, '%s\n', ['PZsparse linear_acc_new', num2str(j),'(0);']);
end
fprintf(fid, '\n');

for i = 1:num_joints
    if i == 1
        switch joint_types{i}
            case 'revolute'
                % (6.45) angular velocity
                w_new = R_t(:,:,i) * w0 + joint_vel(robot_params.q_index(i))*z(:,i); % line 13

                % auxillary angular velocity
                w_aux_new = R_t(:,:,i) * w0_aux + joint_vel_aux(robot_params.q_index(i))*z(:,i); % line 13

                % (6.46) angular acceleration
                wdot_new = R_t(:,:,i) * w0dot ... % line 15
                            + cross(R_t(:,:,i) * w0_aux, joint_vel(robot_params.q_index(i))*z(:,i))...
                            + joint_acc(robot_params.q_index(i))*z(:,i);

                % (6.47) linear acceleration        
                linear_acc_new = R_t(:,:,i)*(linear_acc0 ...
                                    + cross(w0dot, P(:,i)) ... % line 16 (TYPO IN PAPER)
                                    + cross(w0, cross(w0, P(:,i))));      
            case 'fixed'
                % (6.45) angular velocity
                w_new = R_t(:,:,i) * w0; % line 13

                % auxillary angular velocity
                w_aux_new = R_t(:,:,i) * w0_aux; % line 13

                % (6.46) angular acceleration
                wdot_new = R_t(:,:,i) * w0dot; % line 15

                % (6.47) linear acceleration        
                linear_acc_new = R_t(:,:,i)*(linear_acc0 ...
                                    + cross(w0dot, P(:,i)) ... % line 16 (TYPO IN PAPER)
                                    + cross(w0, cross(w0, P(:,i))));  
            otherwise
                error('Joint type not implemented');
        end
                        
    else
        switch joint_types{i}
            case 'revolute'
                % (6.45) angular velocity
                w_new = R_t(:,:,i) * w_old + joint_vel(robot_params.q_index(i))*z(:,i); % line 13

                % auxillary angular velocity
                w_aux_new = R_t(:,:,i) * w_aux_old + joint_vel_aux(robot_params.q_index(i))*z(:,i); % line 14

                % (6.46) angular acceleration
                wdot_new = R_t(:,:,i) * wdot_old ... % line 15
                                    + cross(R_t(:,:,i) * w_aux_old, joint_vel(robot_params.q_index(i))*z(:,i))...
                                    + joint_acc(robot_params.q_index(i))*z(:,i);

                % (6.47) linear acceleration
                linear_acc_new = R_t(:,:,i)*(linear_acc_old ...
                                        + cross(wdot_old, P(:,i)) ... % line 16 (TYPO IN PAPER)
                                        + cross(w_old, cross(w_aux_old,P(:,i))));
            case 'fixed'
                % (6.45) angular velocity
                w_new = R_t(:,:,i) * w_old; % line 13

                % auxillar angular velocity
                w_aux_new = R_t(:,:,i) * w_aux_old; % line 14

                % (6.46) angular acceleration
                wdot_new = R_t(:,:,i) * wdot_old; % line 15

                % (6.47) linear acceleration
                linear_acc_new = R_t(:,:,i)*(linear_acc_old ...
                                        + cross(wdot_old, P(:,i)) ... % line 16 (TYPO IN PAPER)
                                        + cross(w_old, cross(w_aux_old,P(:,i))));
            otherwise
                error('Joint type not implemented');
        end
    end  

    fprintf(fid, '// joint %d\n', i);
    ccode(w_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', 'w_new');
    fprintf(fid, '%s\n', resstr);
    ccode(w_aux_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', 'w_aux_new');
    fprintf(fid, '%s\n', resstr);
    ccode(wdot_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', 'wdot_new');
    fprintf(fid, '%s\n', resstr);
    ccode(linear_acc_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', 'linear_acc_new');
    fprintf(fid, '%s\n', resstr);

    for j = 1:3
        fprintf(fid, '%s\n', ['w', num2str(j),' = w_new',num2str(j),';']);
    end
    for j = 1:3
        fprintf(fid, '%s\n', ['w_aux', num2str(j),' = w_aux_new',num2str(j),';']);
    end
    for j = 1:3
        fprintf(fid, '%s\n', ['wdot', num2str(j),' = wdot_new',num2str(j),';']);
    end
    for j = 1:3
        fprintf(fid, '%s\n', ['linear_acc', num2str(j),' = linear_acc_new',num2str(j),';']);
    end
    fprintf(fid, '\n');

    % (6.48) linear acceleration of CoM auxilliary
    linear_acc_com = linear_acc_old ... % line 23 (modified for standard RNEA)
                    + cross(wdot_old,com(:,i)) ...
                    + cross(w_old,cross(w_aux_old,com(:,i)));                

    % (6.49) calculate forces
    F_new = mass(:,i) * linear_acc_com; % line 27

    ccode(F_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', ['PZsparse F',num2str(i),'_']);
    fprintf(fid, '%s\n', resstr);

    % (6.50) calculate torques
    N_new = I{i}*wdot_old ... % calculated in line 29
             + cross(w_aux_old, (I{i}*w_old));

    ccode(N_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', ['PZsparse N',num2str(i),'_']);
    fprintf(fid, '%s\n', resstr);
end

%% RNEA reverse recursion
F = sym('F', [robot_params.num_q,3], 'real');
F = F';
N = sym('N', [robot_params.num_q,3], 'real');
N = N';

f_old = sym('f_old', [3,1], 'real');
n_old = sym('n_old', [3,1], 'real');

for i = num_joints:-1:1
    if i == num_joints
        % (6.51)
        f_new = F(:,i); % line 28
    
        % (6.52)
        n_new = N(:,i) ...
               + cross(com(:,i), F(:,i));
    else
        % (6.51)
        f_new = R(:,:,i+1) * f_old + F(:,i); % line 28
    
        % (6.52)
        n_new = N(:,i) ...
               + R(:,:,i+1) * n_old ... % line 29
               + cross(com(:,i), F(:,i)) ... % P(:,i) might not be right
               + cross(P(:,i+1), R(:,:,i+1)*f_old); % line 29 (TYPO IN PAPER)
    end

    ccode(f_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', ['PZsparse f',num2str(i),'_'], [num2str(i+1), '_']);
    fprintf(fid, '%s\n', resstr);

    ccode(n_new, 'File', 'temp.txt');
    resstr = handleccodetext('temp.txt', ['PZsparse n',num2str(i),'_'], [num2str(i+1), '_']);
    fprintf(fid, '%s\n', resstr);
end

% calculate joint torquess
for i = 1:num_joints
    switch joint_types{i}
        case 'revolute'
            if z(1,i) == 1
                fprintf(fid, 'u(%d,s_ind) = n%d_1 + damping[%d] * traj->qd_des(%d, s_ind) + armature[%d] * traj->qdda_des(%d, s_ind);\n', robot_params.q_index(i)-1, i, robot_params.q_index(i)-1, robot_params.q_index(i)-1, robot_params.q_index(i)-1, robot_params.q_index(i)-1);
            elseif z(2,i) == 1
                fprintf(fid, 'u(%d,s_ind) = n%d_2 + damping[%d] * traj->qd_des(%d, s_ind) + armature[%d] * traj->qdda_des(%d, s_ind);\n', robot_params.q_index(i)-1, i, robot_params.q_index(i)-1, robot_params.q_index(i)-1, robot_params.q_index(i)-1, robot_params.q_index(i)-1);
            else
                fprintf(fid, 'u(%d,s_ind) = n%d_3 + damping[%d] * traj->qd_des(%d, s_ind) + armature[%d] * traj->qdda_des(%d, s_ind);\n', robot_params.q_index(i)-1, i, robot_params.q_index(i)-1, robot_params.q_index(i)-1, robot_params.q_index(i)-1, robot_params.q_index(i)-1);
            end
        case 'fixed'
            continue;
        otherwise
            error('Joint type not implemented');
    end
end

fclose(fid);

%% post processing
fid = fopen(filename, 'r');
    
resstr = '';
tline = fgetl(fid);
while ischar(tline)
    % replace mass
    for i = 1:robot_params.num_q
        tline = strrep(tline, ['mass',num2str(i)], ['mass_arr(',num2str(i-1),',0)']);
    end

    % replace com
    for i = 1:robot_params.num_q
        for j = 1:3
            tline = strrep(tline, ['com',num2str(j),'_',num2str(i)], ['com[',num2str(i-1),'][',num2str(j-1),']']);
        end
    end

    for i = 1:robot_params.num_q
        tline = strrep(tline, ['I',num2str(i),'_1_1'], ['I_arr(',num2str(i-1),',0)(0,0)']);
        tline = strrep(tline, ['I',num2str(i),'_1_2'], ['I_arr(',num2str(i-1),',0)(0,1)']);
        tline = strrep(tline, ['I',num2str(i),'_1_3'], ['I_arr(',num2str(i-1),',0)(0,2)']);
        tline = strrep(tline, ['I',num2str(i),'_2_1'], ['I_arr(',num2str(i-1),',0)(1,0)']);
        tline = strrep(tline, ['I',num2str(i),'_2_2'], ['I_arr(',num2str(i-1),',0)(1,1)']);
        tline = strrep(tline, ['I',num2str(i),'_2_3'], ['I_arr(',num2str(i-1),',0)(1,2)']);
        tline = strrep(tline, ['I',num2str(i),'_3_1'], ['I_arr(',num2str(i-1),',0)(2,0)']);
        tline = strrep(tline, ['I',num2str(i),'_3_2'], ['I_arr(',num2str(i-1),',0)(2,1)']);
        tline = strrep(tline, ['I',num2str(i),'_3_3'], ['I_arr(',num2str(i-1),',0)(2,2)']);
    end

    resstr = [resstr, tline, newline];

    tline = fgetl(fid);
end

fclose(fid);

fid = fopen(filename, 'w');

fprintf(fid, '%s', resstr);

fclose(fid);

%% helper functions
function resstr = handleccodetext(filename, varName, inputVarName)
    fid = fopen(filename, 'r');
    
    resstr = '';
    tline = fgetl(fid);
    while ischar(tline)
        tline = tline(3:end);

        tline = strrep(tline, 'A0[', varName);
        tline = strrep(tline, '][0]', '');
        
        if nargin < 3
            tline = strrep(tline, '_old', '');
        else
            tline = strrep(tline, '_old', inputVarName);
        end

        if tline(1) ~= 't'
            tline = strrep(tline, '2 = ', '3 = ');
            tline = strrep(tline, '1 = ', '2 = ');
            tline = strrep(tline, '0 = ', '1 = ');
        end

        resstr = [resstr, tline, newline];
        tline = fgetl(fid);
    end

    fclose(fid);
end
