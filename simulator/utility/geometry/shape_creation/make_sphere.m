function varargout = make_sphere(r,c,n)
% S = make_sphere(r,c,n)
% [F,V] = make_sphere(r,c,n)
%
% Make n points on a 3-D sphere of radius r centered at c \in \R^3. These
% points are evenly distributed over the surface of the sphere using
% Vogel's method as outlined in [1,2].
%
% All inputs are optional. The defaults are:
%   r = 1, c = zeros(3,1), n = 100
%
% The number of outputs is optional:
%   S is a 3-by-n (approximately) array of points on the sphere
%   [F,V] is a faces/vertices pair to be used with patch
%
% References:
%   [1] https://bduvenhage.me/geometry/2019/07/31/generating-equidistant-vectors.html
%   [2] http://blog.marmakoide.org/?p=1
%
% Authors: Shreyas Kousik
% Created: 15 Mar 2021
% Updated: 17 Mar 2021

if nargin < 1
    r = 1 ;
end

if nargin < 2
    c = zeros(3,1) ;
end

if nargin < 3
    n = 256 ;
end

% compute golden angle per http://blog.marmakoide.org/?p=1
ga = pi*(3 - sqrt(5)) ;

% create theta angles
T = [0,0,ga * (0:(n-3))] ;

% create z values
Z = [-1,1,linspace(1 - (1/n), (1/n) - 1, n - 2)] ;

% create radii
R = sqrt(1 - Z.^2) ;

% create points at the given cylindrical coordinates
S = [R.*cos(T) ;
    R.*sin(T) ;
    Z] ;

% dilate points by r and shift to c
S = r.*S + repmat(c,1,n) ;

% create output
switch nargout
    case 1
        varargout = {S} ;
    case 2
        S_1 = S(1,:)' ;
        S_2 = S(2,:)' ;
        S_3 = S(3,:)' ;
        F = convhull(S_1,S_2,S_3) ;
        varargout = {F,S'} ;
end