function [q_traj, qd_traj, qdd_traj, t_traj] = get_parameterized_trajectory(q_0, qd_0, k)
    % assuming using ARMTD dynamics...
    t_to_peak = 0:0.005:0.5;
    t_plan = 0.5;
    t_to_stop = 0.505:0.005:1;
    
    for i = 1:length(t_to_peak)
        q_to_peak(:, i) = q_0 + qd_0*t_to_peak(i) + 0.5*k*t_to_peak(i)^2;
        qd_to_peak(:, i) = qd_0 + k*t_to_peak(i);
        qdd_to_peak(:, i) = k;
    end
    
    q_peak = q_to_peak(:, end);
    qd_peak = qd_to_peak(:, end);
    
    braking_acceleration = (0 - qd_peak)/0.5; %final velocity of 0, gets there in 0.5 seconds.
    
    for i = 1:length(t_to_stop)
       q_to_stop(:, i) = q_peak + qd_peak*(t_to_stop(i) - t_plan) + 0.5*braking_acceleration*(t_to_stop(i) - t_plan)^2;
       qd_to_stop(:, i) = qd_peak + braking_acceleration*(t_to_stop(i) - t_plan);
       qdd_to_stop(:, i) = braking_acceleration;
    end
    
    t_traj = [t_to_peak, t_to_stop];
    q_traj = [q_to_peak, q_to_stop];
    qd_traj = [qd_to_peak, qd_to_stop];
    qdd_traj = [qdd_to_peak, qdd_to_stop];
    
    % downsample to get same # of time steps as CORA
    t_traj = t_traj(2:2:end);
    q_traj = q_traj(:, 2:2:end);
    qd_traj = qd_traj(:, 2:2:end);
    qdd_traj = qdd_traj(:, 2:2:end);

end