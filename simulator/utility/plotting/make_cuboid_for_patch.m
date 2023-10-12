function varargout = make_cuboid_for_patch(l,w,h,c)
% [F,V] = make_cuboid_for_patch([l,w,h])
% [F,V] = make_cuboid_for_patch(l,w,h)
% [F,V] = make_cuboid_for_patch(l,w,h,c)
% patch_data_struct = make_cuboid_for_patch(...)
%
% Generate 6 faces and 8 vertices of an axis-aligned cuboid (in 3-D) of
% length l in x, width w in y, and height h in Z. An optional fourth
% argument c defines the center of the cuboid, c \in \R^3/
%
% The output of this function can be passed directly into the patch
% function as faces and vertices, to plot a cuboid in 3-D.
%
% using https://www.mathworks.com/help/matlab/visualize/multifaceted-patches.html
%
% Author: Shreyas Kousik
% Updated: 23 Apr 2020

    if nargin < 4
        c = zeros(3,1) ;
    end

    if nargin == 1
        dims = l ;
        l = dims(1) ;
        w = dims(2) ;
        h = dims(3) ;
    elseif nargin < 1
        l = 1 ;
        w = 1 ;
        h = 1 ;
    end

    % make vertices
    Vx = l.*[0 1 1 0 0 1 1 0]' - l/2 + c(1) ;
    Vy = w.*[0 0 1 1 0 0 1 1]' - w/2 + c(2) ;
    Vz = h.*[0 0 0 0 1 1 1 1]' - h/2 + c(3) ;
    V = [Vx Vy Vz] ;

    % make faces
    F = [1 2 6 5 ; 2 3 7 6 ; 3 4 8 7 ; 4 1 5 8 ; 1 2 3 4 ; 5 6 7 8 ] ;
    
    switch nargout
        case 1
            patch_data_struct.faces = F ;
            patch_data_struct.vertices = V ;
            varargout = {patch_data_struct} ;
        case 2
            varargout = {F,V} ;
        otherwise
            error('This function returns either one or two outputs!')
    end
    
end