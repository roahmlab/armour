function robot = add_ee_link(robot)
    % get robot parameters
    params = get_robot_params(robot);
    
    % create rigid body for the end-effector
    ee_link = robotics.RigidBody('ee_link');
    
    % set mass and inertio to 0
    ee_link.Mass = 0;
    ee_link.Inertia = zeros(1, 6);
    
    % create fixed joint
    ee_link.Joint = robotics.Joint('joint_ee', 'fixed');
    ee_link.Joint.setFixedTransform(eye(4));
    
    % add end-effector to robot
    robot.addBody(ee_link, robot.Bodies{params.active_joints(end)}.Name);
end