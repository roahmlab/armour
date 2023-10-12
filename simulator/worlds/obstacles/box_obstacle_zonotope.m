classdef box_obstacle_zonotope < box_obstacle
    properties
        zono
        Z
        is_base_obstacle = false ;
        creation_buffer = 0 ; % buffer distance when creating random obstacles
        
        buffered_plot_patch_data ;
        buffered_collision_check_patch_data ;
    end
    
    methods
        function O = box_obstacle_zonotope(varargin)
            O@box_obstacle(varargin{:}) ;
            O.make_zono()
            
            O.create_buffered_plot_patch_data ;
            O.create_buffered_collision_check_patch_data ;
        end
        
        function make_zono(O)
            C = O.center ;
            S = O.side_lengths ;
            O.zono = zonotope([C,diag(S./2)]) ;
            O.Z = O.zono.Z;
        end
        
        %% shift center
        function change_center(O,new_center)
            change_center@box_obstacle(O,new_center) ;
            
            % update zonotope
            O.make_zono() ;
        end
        
        %% buffered plot and collision check setup for obstacle generation
        function create_buffered_plot_patch_data(O)
            l = O.side_lengths + 2*O.creation_buffer ;
            c = O.center ;
            
            switch O.dimension
                case 2
                    [F,V] = make_box(l,c) ;
                case 3
                    [F,V] = make_cuboid_for_patch(l(1),l(2),l(3),c) ;
            end
            
            O.buffered_plot_patch_data.faces = F ;
            O.buffered_plot_patch_data.vertices = V ;
        end
        
        function create_buffered_collision_check_patch_data(O)            
            switch O.dimension
                case 2
                    O.buffered_collision_check_patch_data.faces = O.buffered_plot_patch_data.faces ;
                    O.buffered_collision_check_patch_data.vertices = O.buffered_plot_patch_data.vertices ;
                case 3
                    V = O.buffered_plot_patch_data.vertices ;
                    F = convhull(V) ;
                    O.buffered_collision_check_patch_data.faces = F ;
                    O.buffered_collision_check_patch_data.vertices = V ;
            end
        end
        
    end
end