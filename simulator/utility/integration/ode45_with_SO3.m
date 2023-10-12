function [t_out,y_out,R_out] = ode45_with_SO3(dyn,t_span,y_0,R_0,O_idxs,options)
% [t_out,y_out,R_out] = ode45_with_SO3(dyn,t_span,y_0,R_0,O_idxs)
% [t_out,y_out,R_out] = ode45_with_SO3(dyn,t_span,y_0,R_0,O_idxs,options)
%
% This function runs ode45 integration with a rotation matrix as an
% additional state.  The inputs are similar to those of standard MATLAB ode
% solvers, except that the dynamics takes in a rotation matrix, and an
% initial rotation matrix should be passed in as an argument.
%
% INPUTS:
%
%   dyn(t,y,R) is a function that returns the dynamics of the states y
%   
%   t_span = [t_start, t_end] can be an interval or vector of times
%
%   y_0 is an n_states-by-1 vector of initial conditions; n_states must be
%       at least 3, since we must at least have the angular velocity as
%       states of the system to be integrated
%   
%   R_0 is a 3-by-3 initial rotation matrix for the system
%
%   O_idxs is a vector of 3 integers such that y(O_idxs) is the angular
%       velocity, which is needed to compute the dynamics of the rotation
%       matrix R ( as dR/dt = R*skew(O))
%
%   options is the standard ode45 options; this is an optional (hehe)
%       argument to this function
%
% OUTPUTS:
%   t_out is a n_t-by-1 vector of times (if t_span is a vector, then t_out
%       is equal to t_span)
%
%   y_out is an n_t-by-n_states array of states corresponding to each time
%       in t_out
%
%   R_out is a 3-by-3-by-n_t array of rotation matrices, each corresponding
%       to a time in t_out
%
% EXAMPLE USAGE:
%   % angular velocity
%   O_fun = @(t) [sin(t) ; cos(t) ; sin(t)] ;
% 
%   % time
%   t_span = [0,2] ;
%   
%   % initial condition
%   y_0 = O_fun(0) ;
%   R_0 = eye(3) ;
%   O_idxs = 1:3 ;
%
%   % dynamics
%   dyn = @(t,y,R) O_fun(t) ;
%
%   [t_out,y_out,R_out] = ode45_with_SO3(dyn,t_span,y_0,R_0,O_idxs) ;
%
% SEE ALSO: ode45, ode1_with_SO3, ode2_with_SO3
%
% Author: Shreyas Kousik
% Created: 11 Aug 2020
% Updated: -

    % create initial condition
    y_0 = y_0(:) ;
    z_0 = [y_0 ; R_0(:)] ;
    n_states = length(y_0) ;
    
    % create augmented dynamics
    dyn_in = @(t,z) aug_dyn(t,z,dyn,n_states,O_idxs) ;
    
    % call solver
    if nargin > 5
        [t_out,z_out] = ode45(dyn_in, t_span, z_0, options) ;
    else
        [t_out,z_out] = ode45(dyn_in, t_span, z_0) ;
    end
    
    % extract non-rotation and rotation matrix states
    y_out = z_out(:,1:n_states) ;
    R_out = reshape(z_out(:,(n_states+1):end)',3,3,length(t_out)) ;
end

function zd = aug_dyn(t,z,dyn,n_states,O_idxs)
% zd = aug_dyn(t,z,dyn,O_idxs)
%
% Augment the user-provided dynamics with the rotation matrix dynamics.


    % get the non-rotation matrix states, angular velocity state, and
    % rotation matrix state
    y = z(1:n_states) ;
    O = y(O_idxs) ;
    R = reshape(z((end-8):end),3,3) ;
    
    % compute the dynamics for y and R
    dydt = dyn(t,y,R) ;
    dRdt = R*skew(O) ;
    
    % output the augmented dynamics
    zd = [dydt(:) ; dRdt(:)] ;
end

function S = skew(v)
% S = skew(v)
%
% Compute the skew symmetric matrix S from the vector v (this is
% commonly called the "hat map"
    S = [ 0   -v(3)  v(2) ;
         v(3)   0   -v(1) ;
        -v(2)  v(1)   0   ];
end