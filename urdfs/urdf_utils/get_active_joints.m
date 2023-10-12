function joint_idx = get_active_joints(robot)
    for i = 1:robot.NumBodies
        if ~strcmp(robot.Bodies{i}.Joint.Type, 'fixed')
            joint_idx = [joint_idx i];
        end
    end
end