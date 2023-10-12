function shelf_cell = make_shelf_obstacle(center,height,width,depth,...
    N_shelves,min_shelf_height,max_shelf_height,shelf_direction)
% shelf_obs = make_shelf_obstacle(center,height,width,depth,N_shelves,
%                                 min_shelf_height,max_shelf_height,
%                                 shelf_direction)
%
% Creates a cell array of box_obstacle_zonotope objects that represent a
% shelf. The shelf is centered at center \in \R^3. Its height, width, and
% depth depend on the shelf_direction, which can be 1 (the openings on the
% shelf are facing the 'x' direction) or 2 (facing the 'y' direction). The
% shelf has N_shelves shelves, distributed evenly between the min and max
% heights. Each shelf, and the shelf's sides, are 1 cm thick.
%
% Author: Shreyas Kousik
% Created: 23 Dec 2019
% Updated: nah

    % determine the shelf direction
    if nargin < 8
        shelf_direction = 1 ;
    end

    % default shelf thiccness
    thickness = 0.01 ; % 1 cm

    %% figure out which way the shelf is facing
    % 1 for x, 2 for y
    switch shelf_direction
        case 1
            center_1 = center + [0;-width/2;0] ;
            center_2 = center + [0;+width/2;0] ;
            side_side_lengths = [depth thickness height] ;
            shelf_side_lengths = [depth width thickness] ;
        case 2
            center_1 = center + [-width/2;0;0] ;
            center_2 = center + [+width/2;0;0] ;
            side_side_lengths = [thickness depth height] ;
            shelf_side_lengths = [width depth thickness] ;
        otherwise
            error('Please pick 1 or 2 as the shelf direction; 1 is x-facing and 2 is y-facing')
    end

    %% make the sides
    side_1 = {box_obstacle_zonotope('center',center_1,...
                'side_lengths',side_side_lengths,'creation_buffer',0.02)} ;

    side_2 = {box_obstacle_zonotope('center',center_2,...
                'side_lengths',side_side_lengths,'creation_buffer',0.02)} ;

    %% make the shelves
    if nargin < 6
        min_shelf_height = center(3) - height/2 ;
        max_shelf_height = center(3) + height/2 ;
    end

    shelf_heights = linspace(min_shelf_height,...
        max_shelf_height,...
        N_shelves) ;
    shelves = cell(1,N_shelves) ;

    for idx = 1:N_shelves
        shelves{idx} = box_obstacle_zonotope('center',[center(1:2) ; shelf_heights(idx)],...
            'side_lengths',shelf_side_lengths,'creation_buffer',0.02) ;
    end

    %% create output
    shelf_cell = [side_1, side_2, shelves] ;
end