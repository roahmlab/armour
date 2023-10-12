function u = rnea_passivity(q, qd, qd_ref, qdd_ref, robot_params)
    gravity_flag = 1;
    u = rnea(q, qd, qd_ref, qdd_ref, gravity_flag, robot_params);
end