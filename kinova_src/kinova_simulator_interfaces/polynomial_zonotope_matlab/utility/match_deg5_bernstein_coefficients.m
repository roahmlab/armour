function [beta] = match_deg5_bernstein_coefficients(traj_constraints, T)
    % match coefficients to initial position, velocity, acceleration (t=0)
    % and final position, velocity, and acceleration (t=1)
    % assuming a 5th degree bernstein polynomial (minimum degree necessary
    % to satisfy these constraints)
    % also assuming t \in [0, T]

    if nargin < 2
        T = 1;
    end
    
    q0 = traj_constraints{1};
    qd0 = traj_constraints{2};
    qdd0 = traj_constraints{3};
    q1 = traj_constraints{4};
    qd1 = traj_constraints{5};
    qdd1 = traj_constraints{6};
    
    % position constraints:
    beta{1} = q0; % beta_0;
    beta{6} = q1; % beta_5;
    
    % velocity constraints:
    beta{2} = q0 + (T*qd0)/5; % beta_1
    beta{5} = q1 - (T*qd1)/5; % beta_4
    
    % acceleration constraints
    beta{3} = (qdd0*T^2)/20 + (2*qd0*T)/5 + q0; % beta_2
    beta{4} = (qdd1*T^2)/20 - (2*qd1*T)/5 + q1; % beta_3 
    
end