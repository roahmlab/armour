function [t_out,y_out,R_out] = ode2_with_SO3(dyn,t_span,y_0,R_0,O_idxs,dt)
% [t_out,y_out,R_out] = ode2_with_SO3(dyn,t_span,y_0,R_0,O_idxs,dt)
%
% This implements 2nd order RK integration with a rotation matrix state, so
% we're integrating directly on SO(3). The inputs are similar to those of
% standard MATLAB ode solvers, except that the dynamics takes in a rotation
% matrix, and an initial rotation matrix should be given as an argument to
% this function.
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
%       matrix R (as dR/dt = R*skew(O))
%
%   dt is the time discretization to be used (this input is optional, and
%       is ignored if t_span is an array of times
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
%   [t_out,y_out,R_out] = ode2_with_SO3(dyn,t_span,y_0,R_0,O_idxs) ;
%
% SEE ALSO: ode45, ode45_with_SO3, ode1_with_SO3
%
% Author: Shreyas Kousik
% Created: shrug
% Updated: 11 Aug 2020

    %% parse inputs
    if nargin < 5
        % assume last three indices of y are the angular velocity (3-by-1
        % vector that can be mapped to so(3))
        N_states = length(y_0) ;
        O_idxs = (N_states-2):N_states ;
    end
    
    if nargin < 6
        % take ten time steps by default
        dt = (t_span(end) - t_span(1))/10 ;
    end
    
    %% setup
    % setup timing
    if length(t_span) > 2
        t_out = t_span ;
    else
        t_out = t_span(1):dt:t_span(end) ;
        if t_out(end) < t_span(end)
            t_out = [t_out, t_span(end)] ;
        end
    end
    dt_vec = diff(t_out) ;
    N_steps = length(t_out) ;
    
    % allocate output
    y_out = [y_0(:), nan(length(y_0),N_steps-1)] ;
    R_out = cat(3,R_0,nan(3,3,N_steps-1)) ;
    
    %% run integration loop
    for idx = 2:N_steps
        % get time and time step
        t_idx = t_out(idx) ;
        dt_idx = dt_vec(idx-1) ;
        half_dt_idx = dt_idx/2 ;
        
        % get current state and orientation
        y1_idx = y_out(:,idx-1) ;
        R1_idx = R_out(:,:,idx-1) ;
        
    % compute state and rotation matrix dynamics with RK2
        % get first time step
        y1_dot = dyn(t_idx,y1_idx,R1_idx) ;
        k1 = half_dt_idx*y1_dot ;
        O1_idx = y1_idx(O_idxs) ;
        F1_idx = expm_SO3(half_dt_idx.*O1_idx) ;
        R2_idx = F1_idx*R1_idx ;

        % get second time step
        y2_dot = dyn(t_idx + half_dt_idx, y1_idx + k1, R2_idx) ;
        if ~any(isnan(y2_dot))
            k2 = dt_idx*y2_dot ;
            y2_idx = y1_idx + k1 ;
            O2_idx = y2_idx(O_idxs) ;
            F2_idx = expm_SO3(half_dt_idx.*(O1_idx + O2_idx)) ;
        else
            warning('NaNs detected in dynamics! Attempting Euler step.')
            k2 = dt_idx*y1_dot ;
            F2_idx = expm_SO3(dt_idx.*O1_idx) ;
        end

        % compute new state
        y_out(:,idx) = y1_idx + k2 ;
        R_out(:,:,idx) = F2_idx*R1_idx ;
    end
    
    %% fix output
    t_out = t_out(:) ;
    y_out = y_out' ;
end

%% helper functions
function S = skew(O)
    % compute the skew symetric matrix S from angular velocity O
    S = [ 0   -O(3)  O(2) ;
         O(3)   0   -O(1) ;
        -O(2)  O(1)   0   ];
end


function R = expm_SO3(w)
% R = expm_SO3(w)
%
% Computes the vectorized exponential map (at the identity) as defined in:
% http://ethaneade.com/lie.pdf
%
% This version is from: https://github.com/RossHartley/lie.

    theta = vecnorm(w);
    if theta == 0
        R = eye(3);
    else
        w = skew(w) ;
        R = eye(3) + (sin(theta)/theta)*w + ((1-cos(theta))/(theta^2))*w^2;
    end
end