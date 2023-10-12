function [R_out, p_out] = pzfk(R_in, robot_params)
% pass in PZ rotation matrices and robot parameters
% output frame transforms from body frame to world frame for each body

R_out = cell(robot_params.num_joints, 1);
p_out = cell(robot_params.num_joints, 1);

P = robot_params.P;
R0 = matPolyZonotope_ROAHM(eye(3)); % assume base frame is world frame;

for i = 1:robot_params.num_joints
   switch robot_params.joint_types{i}
        case 'revolute'
            R{i, 1} = matPolyZonotope_ROAHM(robot_params.T0(1:3, 1:3, i)) * R_in{robot_params.q_index(i), 1};
        case 'fixed'
            R{i, 1} = matPolyZonotope_ROAHM(robot_params.T0(1:3, 1:3, i));
        otherwise
            error('Joint type not implemented');
   end
   
   if (i == 1)
       p_out{i, 1} = polyZonotope_ROAHM(P(:, i));
       R_out{i, 1} = R0*R{i, 1};
   else
       p_out{i, 1} = p_out{i-1} + R_out{i-1, 1}*P(:, i);
       p_out{i, 1} = reduce(p_out{i, 1}, 'girard', robot_params.zono_order);
       R_out{i, 1} = R_out{i-1, 1}*R{i, 1};
       R_out{i, 1} = reduce(R_out{i, 1}, 'girard', robot_params.zono_order);
   end
end

end

