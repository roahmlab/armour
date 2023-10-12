function [active, fixed] = get_joint_types(robot)
    active = [];
    fixed = [];
    for i = 1:robot.NumBodies
        if strcmp(robot.Bodies{i}.Joint.Type, 'fixed')
            fixed = [fixed i];
        else
            active = [active i];
        end
    end
end