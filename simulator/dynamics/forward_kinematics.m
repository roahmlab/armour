function T = forward_kinematics(q, T0, joint_axes)

    T = eye(4);

    for i =1:length(q)
        dim = find(joint_axes(:,i) ~=0);
        if dim == 1
            R = rx(q(i));
        elseif dim == 2
            R = ry(q(i));
        else
            R = rz(q(i));
        end

        T = T * T0(:,:,i) * [R [0;0;0]; 0 0 0 1];
    end
end