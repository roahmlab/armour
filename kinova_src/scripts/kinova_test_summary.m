%% description
% This script gets the statics of the summaries generated from
% kinova_run_100_worlds
%
% Authors: Bohao Zhang (adapted from Shreyas Kousik code)
% Created 25 October 2022

clear; clc;

use_robust_input = true;

save_file_header = 'trial_' ;
file_location = '../results/random' ;
% file_location = '../results/random_24pi' ;
% file_location = '../results/random_armtd' ;
% file_location = '../results/random_turnoffinputconstraints' ;
% file_location = '../results/random_turnoffinputconstraints_24pi' ;
% file_location = '../results/hard' ;
% file_location = '../results/hard_SL' ;
% file_location = '../results/hard_armtd' ;
% file_location = '../results/hard_armtd_SL' ;
% file_location = '../results/hard_turnoffinputconstraints' ;
addpath(file_location);

summary_files = dir([file_location, '/trial_*']);

collision_check = [];
input_check = [];
ultimate_bound_check = [];
joint_limit_check = [];
goal_check = [];
infeasible_check = [];

for i = 1:length(summary_files)
    data = load(summary_files(i).name);
    if data.summary.collision_check
        collision_check = [collision_check, i];
        continue;
    end
    if data.summary.input_check
        input_check = [input_check, i];
        continue;
    end
    if data.summary.ultimate_bound_check
        ultimate_bound_check = [ultimate_bound_check, i];
        continue;
    end
    if data.summary.joint_limit_check
        joint_limit_check = [joint_limit_check, i];
        continue;
    end
    if data.summary.goal_check
        goal_check = [goal_check, i];
        continue;
    end
    infeasible_check = [infeasible_check, i];
end

fprintf("Test Summary\n");
fprintf("Total Number of Test Trials: %d\n", length(summary_files));
fprintf("Number of Test Trials that collision occurs: %d\n", length(collision_check));
fprintf("Number of Test Trials that exceed torque limits: %d\n", length(input_check));
fprintf("Number of Test Trials that tracking error exceeds ultimate bound: %d\n", length(ultimate_bound_check));
fprintf("Number of Test Trials that exceed joint (position/velocity) limits: %d\n", length(joint_limit_check));
fprintf("Number of Test Trials that reach the goals: %d\n", length(goal_check));
fprintf("Number of Test Trials that do not reach the goals but stop safely: %d\n", length(infeasible_check));
