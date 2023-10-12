function [l,w,h,c] = resize_box_for_world_bounds(l,w,h,c,B)
    % [l,w,h,c] = resize_box_for_world_bounds(l,w,h,c,B)
    %
    % Given a box defined by length/width/height/center and world bounds B
    % given as [xmin xmax ymin ymax zmin zmax], trim the box and recenter
    % it to fit inside the world bounds. If the box lives outside the world
    % bounds, this function throws a warning.
    %
    % See also: bounds_to_box, box_to_bounds
    %
    % Author: Shreyas Kousik
    % Created: ages ago!
    % Updated: 29 Oct 2020
    
    % get the bounds version of the input box
    b = box_to_bounds(l,w,h,c) ;
    
    % adjust the lower bounds
    lo = max([b([1 3 5]) ; B([1 3 5])],[],1) ;
    
    % adjust the upper bounds
    hi = min([b([2 4 6]) ; B([2 4 6])],[],1) ;
    
    % create output
    if any(lo(:) > hi(:))
        l = [] ;
        w = [] ;
        h = [] ;
        c = [] ;
        % warning('Box lies outside of world bounds!')
    else
        [l,w,h,c] = bounds_to_box(reshape([lo;hi],1,6)) ;
    end
end