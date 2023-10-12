function M =  rnea_mass(q, robot_params)
    gravity_flag = 0;
    n = length(q);
    
%     if ~robot_params.use_interval
    if ~strcmp(robot_params.set_type, 'interval')
        M = zeros(n,n);
    else
        M = interval(zeros(n,n), zeros(n,n));
    end
    V = eye(n);

    % calculate mass matrix
    for i = 1:n
        v = V(:,i);
        M(:,i) = rnea(q, zeros(1,n), zeros(1,n), v, gravity_flag, robot_params);
    end
