%% Testing differences between the original urdf and the approximated urdf
% 
clear; clc; close all;

%% import robot
robot_original = importrobot('kinova_with_gripper_dumbbell.urdf');
robot_original.Gravity = [0 0 -9.8];
robot_original.DataFormat = 'column';
robot_approx = importrobot('kinova_with_gripper_dumbbell_approx.urdf');
robot_approx.Gravity = [0 0 -9.8];
robot_approx.DataFormat = 'column';

%% testing
err = [];
for i = 1:10000
    q = randomConfiguration(robot_original);
    % figure(1);
    % show(robot_original, q, 'Visuals', 'on', 'Collisions', 'off');
    % figure(2);
    % show(robot_approx, q, 'Visuals', 'on', 'Collisions', 'off');
    
    mass_original = massMatrix(robot_original, q);
    mass_approx = massMatrix(robot_approx, q);
    error = norm(mass_original - mass_approx)/norm(mass_original);
    disp(error)
    err = [err error];
end
error_mean = mean(err);
disp(error_mean);