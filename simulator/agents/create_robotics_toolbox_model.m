function robot = create_robotics_toolbox_model(agent)
% create_robotics_toolbox_model(agent)
%
% Given a robot_arm_agent as the input, create a robotics toolbox model to
% fill in the agent's corresponding property.

    % initialize robot
    robot = robotics.RigidBodyTree;
    
    % get joint info
    J = agent.joint_locations ;
    J_lims = agent.joint_state_limits ;
    J_axes = agent.joint_axes ;

    if agent.dimension == 2
        J = [J(1:2,:) ; zeros(1,agent.n_links_and_joints) ;
            J(3:4,:) ; zeros(1,agent.n_links_and_joints) ] ;
    end
    

    n = agent.n_links_and_joints ;
    for idx = 1:n
        % make link
        new_link_name = ['link_',num2str(idx)] ;
        new_link = robotics.RigidBody(new_link_name);
        switch idx
            case 0
                new_link.Mass = 0 ;
                new_link.Inertia = zeros(1,6) ;
            otherwise
                new_link.Mass = agent.link_masses(idx) ;
                new_link.Inertia = [diag(agent.link_inertia_matrices{idx})', zeros(1,3)] ;
        end
        
        % set link center of mass location
        new_link.CenterOfMass = J(4:6,idx)' ;

        % make joint
        new_joint_name = ['joint_',num2str(idx)] ;
        new_joint = robotics.Joint(new_joint_name,'revolute') ;
        new_joint.HomePosition = 0 ;
        
        switch idx
            case 0
                new_joint.JointAxis = [0 0 1] ;
                new_joint.PositionLimits = [0 0] ;
            otherwise
                new_joint.JointAxis = J_axes(:,idx)' ;
                new_joint.PositionLimits = J_lims(:,idx)' ;
        end
        

        % set joint location on predecessor link
        switch idx
            case 0
                joint_location = [0 0 0] ;
            case 1
                joint_location = J(1:3,idx)' ;
            otherwise
                joint_location = J(1:3,idx)' - J(4:6,idx-1)' ;
        end
        tform_1 = trvec2tform(joint_location) ;
        setFixedTransform(new_joint,tform_1) ;
        new_link.Joint = new_joint ;

        % add link to robot body
        switch idx
            case 1
                parent_name = 'base' ;
            otherwise
                parent_name = ['link_',num2str(idx-1)] ;
        end
        addBody(robot,new_link,parent_name)
    end

    % add end effector
    ee = robotics.RigidBody('end_effector') ;
    ee.Mass = 0 ;
    ee.Inertia = zeros(1,6) ;
    tform_ee = trvec2tform([agent.link_sizes(1,end), 0, 0]) ;
    setFixedTransform(ee.Joint,tform_ee) ;
    addBody(robot,ee,['link_',num2str(n)])

    % set gravity vector
    robot.Gravity = 9.81.*agent.gravity_direction ;

    % set data format to match the simulator
    robot.DataFormat = 'column' ;

    % update agent
    if nargout == 0
        agent.robotics_toolbox_model = robot ;
    end
end