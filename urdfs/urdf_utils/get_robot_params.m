function [params, robot] = get_robot_params(robot_name)
    robot = importrobot(robot_name);
    robot.DataFormat = 'col';
    robot.Gravity = [0 0 -9.81];

    % get joint types
    [active_joints, fixed_joints] = get_joint_types(robot);
    num_joints = length(active_joints);
     
    %% initialize placeholders
    params = {};
    mass = zeros(1, num_joints);
    com = zeros(3, num_joints);
    joint_axes = zeros(3, num_joints);
    
    % transforms when robot is in home configuration
    H = {};
    
    % position of the origin of frame {i} with respect to frame {i-1}
    P = zeros(3, num_joints + 1);
    
    % spatial inertias
    G = {};
    
    % moments of inertia
    I= {};
    
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
        
        if i == 1
            T0(:,:,i) = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Base.Name);
        else
            T0(:,:,i) = getTransform(robot, q0, robot.Bodies{i}.Name, robot.Bodies{i-1}.Name);
        end
        
        % position of the origin of frame {i} with respect to frame {i-1}
        P(:, i) = H{i}(1:3, 4);
        
        % axis of rotation
        joint_axes(:,i) = robot.Bodies{i}.Joint.JointAxis';
             
        % orientation joint frame{i} in space frame {0}
        Tj(:,:,i) = getTransform(robot, zeros(num_joints,1), robot.Bodies{i}.Name, robot.Base.Name);      
        
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
        
        %~~~~~~~~~~~~~~~~~inertial parameters
        % mass
        mass(i) = robot.Bodies{i}.Mass;

        % center of mass
        com(:, i) = robot.Bodies{i}.CenterOfMass';
        
        % inertia
        I{i} = parallel_axis(robot.Bodies{i}.Inertia, mass(i), com(:,i));  
        
        % spatial inertias
        G{i} = [   I{i}        zeros(3,3); 
                zeros(3,3)   mass(i)*eye(3)];        
                       
    end
    
    params.mass = mass;
    params.com = com;
    params.T0 = T0;
    params.H = H;
    params.P = P;
    params.G = G;
    params.I = I;
    params.M = M;
    params.S = S;
    params.num_bodies = robot.NumBodies;
    params.num_joints = num_joints;
    params.active_joints = active_joints;
    params.fixed_joints = fixed_joints;
    params.joint_axes = joint_axes;
    params.gravity = [0 0 -9.81];
    params.use_interval = false;
end