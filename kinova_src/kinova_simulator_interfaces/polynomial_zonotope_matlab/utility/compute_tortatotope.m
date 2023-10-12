function [interval_torques] = compute_tortatotope(q, qd, qd_a, qdd_a, r, params)

    %% get params
    pz_interval_params = params.pz_interval_params ;
    pz_nominal_params = params.pz_nominal_params;
    pz_diff_params = params.pz_diff_params;

    num_steps = params.num_steps;
    num_joints = params.num_joints;
    use_gravity = params.use_gravity;

    k_P = params.k_P;
    ph_P = params.ph_P;
    k_I = params.k_I;
    ph_I = params.ph_I;
    K_r = params.K_r;

    uniform_bound = params.uniform_bound;
    ult_bound = params.ult_bound;
    err_bound = params.err_bound;
    vel_bound = params.vel_bound;

    %% compute tortatotope
    % initialize placeholders
    for i = 1:num_steps    
        % ith tortatotope
        tau_nom{i,1} = cell(num_joints, 1);
        tau_int{i,1} = cell(num_joints, 1);
        tau_dif{i,1} = cell(num_joints, 1);
        disturbance{i,1} = cell(num_joints, 1);
        v{i,1} = cell(num_joints, 1);
        u{i,1} = cell(num_joints, 1);
    end

    % other placeholders
    rho = zeros(num_joints, num_steps);
    k_bound = zeros(num_steps, 1);
    ph_bound = zeros(num_steps, 1);
    Phi_min = zeros(num_joints, num_steps);
    Phi_max = zeros(num_joints, num_steps);
    d = repmat(interval(0,0), num_joints, num_steps);
    interval_tau_nom = repmat(interval(0,0), num_joints, num_steps);
    interval_tau_int = repmat(interval(0,0), num_joints, num_steps);
    interval_tau_dif = repmat(interval(0,0), num_joints, num_steps);
    interval_disturbance = repmat(interval(0,0), num_joints, num_steps);
    interval_robust = repmat(interval(0,0), num_joints, num_steps);
    interval_torques = repmat(interval(0,0), num_joints, num_steps);

    % compute interval torques
    for i = 1:num_steps    
        % passivity-based tortatotope for nominal params
        tau_nom{i,1} = poly_zonotope_rnea(q{i}, qd{i}, qd_a{i}, qdd_a{i}, use_gravity, pz_nominal_params);

        % passivity-based tortatotope for interval params
        tau_int{i,1} = poly_zonotope_rnea(q{i}, qd{i}, qd_a{i}, qdd_a{i}, use_gravity, pz_interval_params);

        % pass-based tortatotope for diff params (computes the disturbance)
        tau_dif{i,1} = poly_zonotope_rnea(q{i}, qd{i}, qd_a{i}, qdd_a{i}, use_gravity, pz_diff_params);

        % compute disturbance
        for j = 1:num_joints
            % the original way converts disturbance to interval before slicing
            d(j,i) = interval(tau_int{i,1}{j,1} + (-1)*tau_nom{i,1}{j,1});

            % compute pz disturbance
            disturbance{i,1}{j,1} = tau_int{i,1}{j,1} + (-1)*tau_nom{i,1}{j,1};   
        end   
        
        % compute interval mass matrix
        Mtmp = pz_to_int_mass(q{i}, pz_interval_params); % CHANGE THIS
        
        M{i} = Mtmp - Mtmp; % CHANGE THIS
        
        % integral action bounds 
        k_bound(i,1) = k_P ;
        ph_bound(i,1) = ph_P ;

        % bounds on max disturbance
        Phi_min(:,i) = abs(infimum(d(:,i)));
        Phi_max(:,i) = abs(supremum(d(:,i)));
        rho(:,i) = max(Phi_min(:,i), Phi_max(:,i));

        for j = 1:num_joints
            % complete tortatotope    
%             v{i,1}{j,1} = (k_bound(i,1)*rho(i,1) + ph_bound(i,1)) * r{i}{j,1};
%             v{i,1}{j,1} = (k_P * norm(rho(:,i)) + ph_P)  * r{i}{j,1};
            v{i,1}{j,1} = 23.875 * norm(M{i}, 'fro') + norm(rho(:,i)) ;
            u{i,1}{j,1} = v{i}{j,1} + tau_nom{i}{j,1};

            % over-approximate terms with intervals
            c(j,i) = interval(tau_nom{i}{j});
            interval_tau_int(j,i) = interval(tau_int{i}{j});
            interval_tau_dif(j,i) = interval(tau_dif{i}{j});
            interval_disturbance(j,i) = interval(disturbance{i}{j});
            interval_robust(j,i) = interval(v{i}{j});
            interval_torques(j,i) = interval(u{i}{j});

            % partition fully sliceable generators
    %         [C{i,1}{j,1}, Gk{i,1}{j,1}, buffer{i,1}{j,1}, E{i,1}{j,1}, ~] = partition_generators(u{i,1}{j,1}, num_joints, id_names);      
        end
    end

end

function M = pz_mass(q, params)
    gravity_flag = 0;
    n = length(q);
   
    % identity poly zonotope
    E = pz_eye(n);
    
    % calculate mass matrix
    for i = 1:n
        % get unit and zero vectors
        e = E{i};
        z = E{i};
        z{i,1} = 0*z{i,1};
        
        % compute ith column
        M{i} = poly_zonotope_rnea(q, z, z, e, gravity_flag, params);
    end
end

function E = pz_eye(n)
    for i = 1:n
        for j = 1:n
            if i == j
                col{j,1} = polyZonotope_ROAHM(1, 0, [], 1, 0);
            else
                col{j,1} = polyZonotope_ROAHM(0, 0, [], 1, 0);
            end
        end
        E{1,i} = col;
    end
end

function M = pz_to_int_mass(q, params)
    gravity_flag = 0;
    n = length(q);
   
    M = interval(zeros(n,n), zeros(n,n));
        
    % identity poly zonotope
    E = pz_eye(n);
    
    % calculate mass matrix
    for i = 1:n
        % get unit and zero vectors
        e = E{i};
        z = E{i};
        z{i,1} = 0*z{i,1};
        
        % compute ith column
        m = poly_zonotope_rnea(q, z, z, e, gravity_flag, params);
        
        % convert to interval
        for j = 1:n
            M(j,i) = interval(m{j,1});
        end
    end
end