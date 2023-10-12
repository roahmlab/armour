classdef arm_manual_waypoint_HLP < manual_waypoint_HLP
    properties
        arm_plot_data_function = [] ;
    end
    
    methods
        function HLP = arm_manual_waypoint_HLP(varargin)
            HLP@manual_waypoint_HLP(varargin{:}) ;
        end
        
        function setup(HLP,agent_info,world_info)
            setup@manual_waypoint_HLP(HLP,agent_info,world_info) ;
            
            % set up for plops
            HLP.arm_plot_data_function = agent_info.get_collision_check_volume ;
        end
        
        function plot(HLP)
            % plot current waypoint
            wp = HLP.current_waypoint ;
            
            if ~isempty(wp)
                D = HLP.arm_plot_data_function(wp) ;
                F = D.faces ;
                V = D.vertices ;
                
                if check_if_plot_is_available(HLP,'current_waypoint')
                    HLP.plot_data.current_waypoint.faces = F ;
                    HLP.plot_data.current_waypoint.vertices = V ;
                else
                    d = patch('faces',F,'vertices',V,...
                        'facealpha',0.1,'facecolor','g',...
                        'edgealpha',0.5,'edgecolor','g') ;
                    HLP.plot_data.current_waypoint = d ;
                end
            end
        end
    end
end