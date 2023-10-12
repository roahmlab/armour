classdef robot_arm_2D_2link_2DOF_thin < robot_arm_agent
    methods
        function A = robot_arm_2D_2link_2DOF_thin(varargin)
            n_links_and_joints = 2 ;
            
            dimension = 2 ;
            
            link_shapes = {'oval','oval'} ;
            
            link_sizes = [1.000, 1.000 ;
                          0.001, 0.001] ;

            joint_locations = [+0.000 +0.500 ;
                               +0.001 +0.000 ;
                               -0.500 -0.500 ;
                                0.000 +0.000 ] ;
            
            joint_state_limits = [-Inf, -Inf ;
                                  +Inf, +Inf] ;
            
            joint_speed_limits = [-Inf, -Inf ;
                                  +Inf, +Inf] ;
            
            joint_input_limits = [-Inf, -Inf ;
                                      +Inf, +Inf ] ;
                                  
            n_states = 2*n_links_and_joints ;
            
            joint_state_indices = 1:2:n_states ;
            
            joint_speed_indices = 2:2:n_states ;
            
            joint_types = {'revolute','revolute','revolute','revolute', 'revolute', 'revolute'} ;
            
            joint_axes = [0 0;
                          0 0;
                          1 1] ; % axes are in preceding link's frame
                      
            kinematic_chain = [0 1;
                               1 2] ;
                           
            gravity_direction = [0;-1;0] ;
            
            set_view_when_animating = false ;
            
            buffer_dist = 0;
            
            A@robot_arm_agent('dimension',dimension,...
                'n_links_and_joints',n_links_and_joints,...
                'n_inputs',n_links_and_joints,...
                'n_states',n_states,...
                'link_shapes',link_shapes,'link_sizes',link_sizes,...
                'joint_state_indices',joint_state_indices,...
                'joint_speed_indices',joint_speed_indices,...
                'joint_types',joint_types,...
                'joint_axes',joint_axes,...
                'joint_locations',joint_locations,...
                'joint_state_limits',joint_state_limits,...
                'joint_speed_limits',joint_speed_limits,...
                'joint_input_limits',joint_input_limits,...
                'kinematic_chain',kinematic_chain,...
                'gravity_direction',gravity_direction,...
                'animation_set_view_flag',set_view_when_animating,...
                'buffer_dist', buffer_dist,...
                'animation_time_discretization', 0.01,...
                varargin{:}) ;
        end
    end
end