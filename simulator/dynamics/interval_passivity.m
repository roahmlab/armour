function u = interval_passivity(q, qd, qd_ref, qdd_ref, robot_params)
    gravity_flag = 1;
    u = interval_rnea(q, qd, qd_ref, qdd_ref, gravity_flag, robot_params);
end