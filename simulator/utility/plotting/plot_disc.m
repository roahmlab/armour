function h = plot_disc(r,p,varargin)
% plot_disc(r,p)
% h = plot_disc(r,p,options)
%
% Plot a disc of radius r centered at p in \R^2. If p is omitted, the disc
% is centered at zero.
%
% This function creates and returns a patch handle. One can pass in the
% regular arguments to the patch function as well.
%
% For now, since I'm lazy, this function just creates a 100-sided polygon
% instead of actually creating a disc.
%
% Authors: Shreyas Kousik
% Created: 7 Jan 2021
% Updated: nah

    % get circle vertices and (single) face
    V = make_circle(r,100,p)' ;
    F = [1:100, 1] ;

    % create patch object
    h = patch('faces',F,'vertices',V,varargin{:}) ;
    
    if nargout < 1
        clear h
    end
end