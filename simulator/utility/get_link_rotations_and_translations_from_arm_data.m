function [R,T,J] = get_link_rotations_and_translations_from_arm_data(j_vals,j_axes,j_locs)
% [R,T,J] = get_link_rotations_and_translations_from_arm_data(j_vals,j_axes,j_locs)
%
% Given joint values, joint axes, and joint locations on each link, return
% the rotation matrices, translation vectors, and joint locations for the
% given configuration (i.e., joint values)
%
% Note, this only works in 3-D right now. The joint axes can be gotten from
% a robot_arm_agent's joint_axes property, and similarly joint locations
% can be gotten from the joint_locations property. The input j_vals should
% a 1-by-n_links_and_joints vector where n_links_and_joints is a property
% of the robot_arm_agent.
%
% Author: Shreyas Kousik
% Created: who knows
% Updated: 29 Feb 2020

    %% setup
    n = length(j_vals) ;
    d = 3 ;
    kinematic_chain = [0:(n-1) ; 1:n] ;

    % allocate cell arrays for the rotations and translations
    R = mat2cell(repmat(eye(d),1,n),d,d*ones(1,n)) ;
    T = mat2cell(repmat(zeros(d,1),1,n),d,ones(1,n)) ;

    % allocate array for the joint locations
    J = nan(d,n) ;
    
    if isa(j_vals,'sym')
        J = sym(J) ;
    end
    
    %% precompute rotation matrices between joint axes
    R_axes = cell(1,n) ;    
    R_axes{1} = eye(3) ; % rotation matrix of first joint axis relative to baselink
    
    for idx = 2:n
        a_1 = j_axes(:,idx-1) ;
        a_2 = j_axes(:,idx) ;
        
        v = cross(a_1,a_2) ;
        s = vecnorm(v) ;
        c = a_1'*a_2 ;
        
        R_axes{idx} = eye(3) + skew(v) + (skew(v)^2).*(1-c)./(s^2) ;
    end

    %% get outputs
    % move through the kinematic chain and get the rotations and
    % translation of each link
    for idx = 1:n
        k_idx = kinematic_chain(:,idx) ;
        p_idx = k_idx(1) ;
        s_idx = k_idx(2) ;

        % get the rotation and translation of the predecessor and
        % successor links; we assume the baselink is always rotated
        % to an angle of 0 and with a translation of 0
        if p_idx == 0
            R_pred = eye(d) ;
            T_pred = zeros(d,1) ;
        else
            R_pred = R{p_idx} ;
            T_pred = T{p_idx} ;
        end

        % get the value and location of the current joint
        j_idx = j_vals(idx) ;
        j_loc = j_locs(:,idx) ;

        % compute link rotation matrix of current link -- this works fine
        % for a numerical representation of the rotation matrix, but things
        % blow up really badly if you plug symbolic variables in here :/
        % -- note left by shreyas 1 oct 2020 --
        axis_pred = R_pred*j_axes(:,idx) ;
        K = skew(axis_pred) ;
        R_succ = eye(3) + sin(j_idx)*K + (1 - cos(j_idx))*(K^2) ;
        R_succ = R_succ*R_pred ;

        % create translation
        T_succ = T_pred + R_pred*j_loc(1:d) - R_succ*j_loc(d+1:end) ;

        % fill in rotation and translation cells
        R{s_idx} = R_succ ;
        T{s_idx} = T_succ ;

        % fill in the joint location
        j_loc_local = j_locs((d+1):end,idx) ;
        J(:,idx) = -R_succ*j_loc_local + T_succ ;
    end
end
