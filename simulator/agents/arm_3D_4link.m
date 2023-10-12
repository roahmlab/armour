classdef arm_3D_4link < robot_arm_agent
    methods
        function A = arm_3D_4link(varargin)
            dimension = 3 ;
            
            link_shapes = {'cuboid', 'cylinder','cuboid','cylinder','cuboid','cylinder','cuboid','cylinder'} ;
            
            link_sizes = [0.01, 0.25, 0.01, 0.25, 0.01, 0.25, 0.01, 0.25;
                0.01, 0.10, 0.01, 0.10, 0.01, 0.10, 0.01, 0.10 ;
                0.01, 0.10, 0.01, 0.10, 0.01, 0.10, 0.01, 0.10] ;
            
            joint_axes = [0 0 1 0 1 0 1 0 ;
                          0 1 0 1 0 1 0 1 ;
                          1 0 0 0 0 0 0 0 ] ;
            
            joint_locations = [0    0      0.25/2 0      0.25/2  0       0.25/2  0 ;
                               0    0      0      0      0       0       0       0 ;
                               0.05 0      0      0      0       0       0       0 ;
                               0   -0.25/2 0     -0.25/2 0      -0.25/2  0      -0.25/2 ;
                               0    0      0      0      0       0       0       0 ;
                               0    0      0      0      0       0       0       0 ] ;
             
            joint_state_limits = [-inf -pi/2 -inf -pi/2 -inf -pi/2 -inf -pi/2;
                                   inf  pi/2  inf  pi/2  inf  pi/2  inf  pi/2];
            
            joint_speed_limits = 2*pi.*repmat([-1;1],1,length(link_shapes)) ;
            
            joint_input_limits = 3.*joint_speed_limits;
            
            A@robot_arm_agent('dimension',dimension,...
                'link_shapes',link_shapes,'link_sizes',link_sizes,...
                'joint_axes',joint_axes,...
                'joint_locations',joint_locations,...
                'joint_state_limits',joint_state_limits,...
                'joint_speed_limits',joint_speed_limits,...
                'joint_input_limits',joint_input_limits,...
                varargin{:}) ;
        end
    end
end