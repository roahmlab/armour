function I = parallel_axis(inertia_vec, mass, com)
    %% PARALLEL AXIS
    % Inputs
    %   inertia_vec - the inertia tensor obtained from RigidBody object
    %   mass - mass of the link
    %   com - the center of mass relative to the joint frame
    %
    % Note that RigidBody expresses the inertia tensor with respect to 
    % the joint frame. This function uses the parallel-axis theorem to 
    % compute the inertia tensor in the link body frame whose origin is
    % located at the link center of mass. This function returns the inertia
    % values that are in the URDF.
    
    % We are computing the following:
    % I_{xx}^{C} = I_{xx}^{A} - m (y_c^2 + z_c^2);
    % I_{yy}^{C} = I_{yy}^{A} - m (x_c^2 + z_c^2);
    % I_{zz}^{C} = I_{zz}^{A} - m (x_c^2 + y_c^2);
    %
    % I_{xy}^{C} = I_{xy}^{A} + m x_c y_c;
    % I_{xz}^{C} = I_{xz}^{A} + m x_c z_c;
    % I_{yz}^{C} = I_{yz}^{A} + m y_c z_c;
    %% compute new inertia tensor
    % mass
    m = mass;

    % center of mass coordinates
    x = com(1);
    y = com(2);
    z = com(3);
    
    % inertia vec from matlab rigidBody object
    Ixx = inertia_vec(1);
    Iyy = inertia_vec(2);
    Izz = inertia_vec(3);
    Iyz = inertia_vec(4);
    Ixz = inertia_vec(5);
    Ixy = inertia_vec(6);
    
    % parallel axis theorem (see Craig 6.25)
    Ixx = Ixx - m * (y^2 + z^2); 
    Iyy = Iyy - m * (x^2 + z^2);
    Izz = Izz - m * (x^2 + y^2);
    
    Ixy = Ixy + m * x * y;
    Ixz = Ixz + m * x * z;
    Iyz = Iyz + m * y * z;
    
    Iyx = Ixy;
    Izx = Ixz;
    Izy = Iyz;
    
    I = [Ixx Ixy Ixz;
         Iyx Iyy Iyz;
         Izx Izy Izz];

end