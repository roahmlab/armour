function [true_params, nominal_params, interval_params, lo_robot, hi_robot] = get_interval_params(robot_name, use_random)
    %% nominal params
    nominal_params = get_robot_params(robot_name);
    
    %% interval params
    [interval_params, lo_robot] = get_robot_params(robot_name);
    [~, hi_robot] = get_robot_params(robot_name);
    
    num_joints = interval_params.num_joints;
    
    lo_mass_scale = 0.97;
    hi_mass_scale = 1.03;
    
    lo_com_scale = 1;
    hi_com_scale = 1;
    
    lo_mass = zeros(1, num_joints);
    hi_mass = zeros(1, num_joints);
  
    lo_com = zeros(3, num_joints);
    hi_com = zeros(3, num_joints);
    
    I = {};
    
    % scale mass, inertia, and com
    for i = 1:nominal_params.num_joints
        % scale link masses
        lo_robot.Bodies{i}.Mass = lo_mass_scale * lo_robot.Bodies{i}.Mass;
        hi_robot.Bodies{i}.Mass = hi_mass_scale * hi_robot.Bodies{i}.Mass;
        lo_mass(i) = lo_robot.Bodies{i}.Mass;
        hi_mass(i) = hi_robot.Bodies{i}.Mass;        
        
        % scale center of mass
        lo_robot.Bodies{i}.CenterOfMass = lo_com_scale * lo_robot.Bodies{i}.CenterOfMass;
        hi_robot.Bodies{i}.CenterOfMass = hi_com_scale * hi_robot.Bodies{i}.CenterOfMass;
        lo_com(:, i) = lo_robot.Bodies{i}.CenterOfMass';
        hi_com(:, i) = hi_robot.Bodies{i}.CenterOfMass';  
        
        % scale inertia vec
        lo_robot.Bodies{i}.Inertia = lo_mass_scale * lo_robot.Bodies{i}.Inertia;
        hi_robot.Bodies{i}.Inertia = hi_mass_scale * hi_robot.Bodies{i}.Inertia;  
        lo_inertia_vec = lo_robot.Bodies{i}.Inertia;
        hi_inertia_vec = hi_robot.Bodies{i}.Inertia;
        
        % create inertia matrix about com 
        lo_I = parallel_axis(lo_inertia_vec,...
                             lo_mass(i),... 
                             lo_com(:, i));
                         
        hi_I = parallel_axis(hi_inertia_vec,...
                             hi_mass(i),... 
                             hi_com(:, i));
                         
        % swap high and low values
        for j = 1:3
            for k = 1:3
                a = min(lo_I(j,k), hi_I(j,k));
                b = max(lo_I(j,k), hi_I(j,k));
                
                lo_I(j,k) = a;
                hi_I(j,k) = b;
            end
        end
        
        % create interval inertia        
        I{i} = interval(lo_I, hi_I); 
          
        a = min(lo_com(:,i), hi_com(:,i));
        b = max(lo_com(:,i), hi_com(:,i));
        
        lo_com(:,i) = a;
        hi_com(:,i) = b;
        
        % interval spatial inertia
        G{i} = [I{i}, interval(zeros(3, 3));....
                interval(zeros(3, 3)), interval(lo_mass(i), hi_mass(i))*eye(3)];
        
    end
     
    % save interval parameters
    interval_params.mass = interval(lo_mass, hi_mass);
    interval_params.com = interval(lo_com, hi_com);
    interval_params.I = I;
    interval_params.G = G;
    interval_params.use_interval = true;
    
    %% true params 
%     true_mass_scale = 1.01;
%     true_com_scale = 0.99;
    
    [true_params, true_robot] = get_robot_params(robot_name);
    for i = 1:num_joints
        if use_random
            true_mass_scale = lo_mass_scale + (hi_mass_scale-lo_mass_scale) .* rand(1,1);
%             true_com_scale = lo_com_scale + (hi_com_scale-lo_com_scale) .* rand(1,1);
        else
            true_mass_scale = 1.01;
%             true_com_scale = 0.99;
        end
        
        true_params.mass(i) = true_mass_scale * true_params.mass(i);
        true_params.I{i} = true_mass_scale * true_params.I{i};
%         true_params.com(:, i) = true_com_scale * ;
        
%         % mass
%         true_robot.Bodies{i}.Mass = true_mass_scale * true_robot.Bodies{i}.Mass;
% 
%         % center of mass
%         true_robot.Bodies{i}.CenterOfMass = true_com_scale * true_robot.Bodies{i}.CenterOfMass;
% 
%         % inertia
%         true_robot.Bodies{i}.Inertia = true_mass_scale * true_robot.Bodies{i}.Inertia;
    end
    
%     true_params = get_robot_params(true_robot);
    
end