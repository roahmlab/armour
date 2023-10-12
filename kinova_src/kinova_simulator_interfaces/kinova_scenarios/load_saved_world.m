function [start, goal, obstacles] = load_saved_world(filename)
% function for loading a presaved world

M = csvread(filename);
start = M(1, :)';
goal = M(2, :)';

obstacles = {};
[rows_M, ~] = size(M);

for i = 4:rows_M
    obstacles{i-3} = box_obstacle_zonotope('center', M(i, 1:3)', 'side_lengths', M(i, 4:6)');
end


end