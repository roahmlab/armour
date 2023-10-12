function X = make_grid_2D(bounds,n_1,n_2)
% X = make_grid_2D(bounds,n_1,n_2)
%
% Make a 2-D grid of points in the bounds [x1lo x1hi x2lo x2hi] with n_1
% samples along the x_1 direction and n_2 samples along the x_2 direction.
% The output X is a 2-by-(n_1*n_2) array of points.
%
% See also: make_grid
%
% Authors: Shreyas Kousik
% Created: 2 Feb 2021 (tbh I made this function, like, eons ago)
% Updated: 9 Feb 2021
    if nargin < 2
        n_1 = 100 ;
    end
    
    if nargin < 3
        n_2 = n_1 ;
    end
    
    if nargin < 1
        bounds = [-1 1 -1 1] ;
    end
    
    x_1 = linspace(bounds(1),bounds(2),n_1) ;
    x_2 = linspace(bounds(1),bounds(2),n_2) ;
    
    [X1,X2] = meshgrid(x_1,x_2) ;
    X = [X1(:), X2(:)]' ;
end