function C = rnea_coriolis(q, qd, robot_params)
    gravity_flag = 0;
    
    n = length(q);
    V = eye(n);
    
%     if ~robot_params.use_interval
    if ~strcmp(robot_params.set_type, 'interval')
        C = zeros(n,n);
    else
        C = interval(zeros(n,n), zeros(n,n));
    end

    % calculate Coriolis matrix
    for i = 1:n
        v = V(:,i);
        C(:,i) = rnea(q, qd, v, zeros(1,n), gravity_flag, robot_params);
    end
end