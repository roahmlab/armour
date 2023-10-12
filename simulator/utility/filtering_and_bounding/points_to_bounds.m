function B = points_to_bounds(P)
% B = points_to_bounds(P)
%
% Given points P as a 2-by-N (resp. 3-by-N) array of points in \R^2 (resp.
% \R^3), return a bounds vector
%   B = [xmin xmax ymin ymax (resp., zmin zmax)]
%
% Author: Shreyas Kousik
% Created: 23 Apr 2020
% Updated: not yet

    B = [min(P,[],2), max(P,[],2)]' ;
    B = B(:)' ;
end