function c = rnea_coriolis_vec(q, qd, robot_params)
    gravity_flag = 0;
    n = length(q);
    c = rnea(q, qd, qd, zeros(1,n), gravity_flag, robot_params);
end