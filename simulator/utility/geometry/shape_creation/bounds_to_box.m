function varargout = bounds_to_box(B)
% varargout = bounds_to_box(B)
%
% Given bounds [xmin xmax ymin ymax zmin zmax], make a box given by
% length/width/height/center [l,w,h,c]
% 
% Usage:
%   [l,w,c] = bounds_to_box([xmin xmax ymin ymax])
%   [l,w,h,c] = bounds_to_box([xmin xmax ymin ymax zmin zmax])
%   V = bounds_to_box([xmin xmax ymin ymax])
%
% V is either a 2-by-5 of 3-by-8 array of points containing the vertices of
% the box
%
% Author: Shreyas Kousik
% Created: Mar 2019
% Updated: 23 Apr 2020
    switch length(B)
        case 4
            B = reshape(B(:),2,2)' ;
            c = mean(B,2) ;
            l = B(1,2) - B(1,1) ;
            w = B(2,2) - B(2,1) ;
            if nargout == 1
                varargout = {make_box([l w],c)} ;
            elseif nargout == 3
                varargout = {l,w,c} ;
            else
                error('Incorrect number of output arguments!')
            end
        case 6
            B = reshape(B(:),2,3)' ;
            c = mean(B,2) ;
            l = B(1,2) - B(1,1) ;
            w = B(2,2) - B(2,1) ;
            h = B(3,2) - B(3,1) ;
            
            if nargout == 1
                % make vertices
                Vx = l.*[0 1 1 0 0 1 1 0]' - l/2 + c(1) ;
                Vy = w.*[0 0 1 1 0 0 1 1]' - w/2 + c(2) ;
                Vz = h.*[0 0 0 0 1 1 1 1]' - h/2 + c(3) ;
                varargout = {[Vx Vy Vz]'} ;
            elseif nargout == 4
                varargout = {l,w,h,c} ;
            else
                error('Incorrect number of output arguments!')
            end
        otherwise
            error(['The input should be either [xmin xmax ymin ymax] ',...
                'or [xmin xmax ymin ymax zmin zmax].'])
    end
end