function u = rnea_invdyn(q, qd, qdd, robot_params)
    gravity_flag = 1;
    u = rnea(q, qd, qd, qdd, gravity_flag, robot_params);
end