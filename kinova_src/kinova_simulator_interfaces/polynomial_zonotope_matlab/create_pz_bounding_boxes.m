function [link_poly_zonotopes, link_sizes, mesh] = create_pz_bounding_boxes(robot)
% takes in a robot
% creates overapproximating PZ bounding box for each body of the robot.

for i = 1:robot.NumBodies
    bounds = zeros(3, 2);

    if ~isempty(robot.Bodies{i}.Visuals)
        stl_file = robot.Bodies{i}.Visuals{1};
        stl_file = extractAfter(stl_file, 'Mesh:');
        mesh{i, 1} = stlread(stl_file);
        
        for j = 1:3
            bounds(j, 1) = min(mesh{i}.Points(:, j));
            bounds(j, 2) = max(mesh{i}.Points(:, j));
        end
    else
        % 10 cm cube if STL file not available
        bounds = [-0.05, 0.05; 
                  -0.05, 0.05; 
                  -0.05, 0.05];
    end
    
    c = (bounds(:, 1) + bounds(:, 2))./2;
    g = (bounds(:, 2) - bounds(:, 1))./2;
    
    link_poly_zonotopes{i, 1} = polyZonotope_ROAHM(c, [], [g(1), 0, 0;...
                                                           0, g(2), 0;...
                                                           0, 0, g(3)]);
                                                       
    link_sizes(:, i) = [2*g(1); 2*g(2); 2*g(3)];
end

end

