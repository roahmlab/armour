function [params, robot] = get_kinematic_params(robot)
    robot.DataFormat = 'col';
    robot.Gravity = [0 0 -9.81];

    %% get joint types
%     [active_joints, fixed_joints] = get_joint_types(robot);
%     num_joints = length(active_joints);
    num_joints = length(robot.Bodies);
    joint_types = cell(num_joints, 1);
    q_index = nan(num_joints, 1);
    q_count = 0;
    for i = 1:robot.NumBodies
       joint_types{i, 1} = robot.Bodies{i}.Joint.Type;
       if ~strcmp(robot.Bodies{i}.Joint.Type, 'fixed')
           % get corresponding configuration for this index.
           q_count = q_count + 1;
           q_index(i, 1) = q_count;
       end
    end
    num_q = q_count;
     
    %% initialize placeholders
    params = {};
    joint_axes = zeros(3, num_joints);
    
    % transforms when robot is in home configuration
    H = {};
    
    % position of the origin of frame {i} with respect to frame {i-1}
    P = zeros(3, num_joints + 1);
    
    % orientation of com frame {i} in com frame {i-1}
    M = repmat(eye(4), [1,1,num_joints+1]);
    
    % screw axes in world frame
    S = zeros(6, num_joints);
    
    % orientation joint frame{i} in space frame {0}
    Tj = zeros(4, 4, num_joints);
    
    % orientation of com frame {i} in space frame {0}
    Tc = zeros(4,4,num_joints);
    
    %% get/compute params
    
    q0 = homeConfiguration(robot);
    T0 = zeros(4,4,num_joints);
    
    for i = 1:num_joints
        %~~~~~~~~~~~~~~~~~transforms 
        % calculate frame-to-frame transforms
        if i == 1
            H{i} = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Base.Name);      
        else
            H{i} = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Bodies{i-1}.Name);
        end   
        
        T0(:, :, i) = H{i};
        
        % position of the origin of frame {i} with respect to frame {i-1}
        P(:, i) = H{i}(1:3, 4);
        
        % axis of rotation
        joint_axes(:,i) = robot.Bodies{i}.Joint.JointAxis';
             
        % orientation joint frame{i} in space frame {0}
        Tj(:,:,i) = getTransform(robot, zeros(num_q,1), robot.Bodies{i}.Name, robot.Base.Name);      
        
        % orientation of com frame {i} in space frame {0}
        K = eye(4);
        K(1:3,4) = robot.Bodies{i}.CenterOfMass' ;
        Tc(:,:,i) = Tj(:,:,i) * K;

        % orientation of com frame {i} in com frame {i-1}
        if i == 1
           M(:,:,i) = Tc(:,:,i);
        else
           M(:,:,i) = TransInv(Tc(:,:,i-1)) * Tc(:,:,i);
        end

        % screw axes
        w = Tj(1:3,1:3,i) * robot.Bodies{i}.Joint.JointAxis' ; 
        v = cross(-w, Tj(1:3, 4, i));
        S(:,i) = [w;v];              
    end
    
    params.T0 = T0;
    params.H = H;
    params.P = P;
    params.M = M;
    params.S = S;
    params.num_bodies = robot.NumBodies;
    params.num_joints = num_joints;
    params.num_q = num_q;
    params.joint_types = joint_types;
    params.q_index = q_index;
    params.joint_axes = joint_axes;
    params.gravity = [0 0 -9.81];

end