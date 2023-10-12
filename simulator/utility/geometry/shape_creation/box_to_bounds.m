function B = box_to_bounds(l,w,h,c)
% B = box_to_bounds(l,w,h,c)
% B = box_to_bounds([l,w,h],c)
%
% Given a box length/width/height/center, make bounds as a vector
%   [xmin xmax ymin ymax zmin zmax]
%
% Note this doesn't work in 2-D!
%
% Author: Shreyas Kousik
% Created: shrug
% Updated: 23 Apr 2020

    if nargin == 2
        c = w ;
        w = l(2) ; 
        h = l(3) ;
        l = l(1) ;
    end

    xmin = c(1) - l/2 ;
    xmax = c(1) + l/2 ;
    ymin = c(2) - w/2 ;
    ymax = c(2) + w/2 ;
    zmin = c(3) - h/2 ;
    zmax = c(3) + h/2 ;
    
   B = [xmin xmax ymin ymax zmin zmax] ;
end