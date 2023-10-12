initialize_script_path = matlab.desktop.editor.getActiveFilename;
cd(initialize_script_path(1:end-12));

%% Initialize controller
fprintf("Start compiling mex robust controller\n\n");

% default robot to test
kinova_model_filename = 'kinova_without_gripper';

robot_model_file_path_filename = [pwd, '/kinova_simulator_interfaces/kinova_robust_controllers_mex/robot_model_file_path.hpp'];
fid = fopen(robot_model_file_path_filename, 'w');
fprintf(fid, 'const char RobotFilePath[] = "%s/kinova_simulator_interfaces/kinova_robust_controllers_mex/%s.txt";', pwd, kinova_model_filename);
fclose(fid);

cd kinova_simulator_interfaces/kinova_robust_controllers_mex/

try
    compile;
catch
    error('Error when compiling mex robust controller code!');
end

fprintf("Successfully compiled mex robust controller\n\n");

cd ../../

%% Initialize real-time armour planner
fprintf("\nStart compiling ARMOUR\n\n");

cd kinova_simulator_interfaces/kinova_planner_realtime

armour_buffer_path = [pwd, '/BufferPath.h'];
fid = fopen(armour_buffer_path, 'w');
fprintf(fid, '#include <string>\n');
fprintf(fid, 'const std::string pathname = "%s/buffer/";', pwd);
fclose(fid);

terminal_output = system('./compile.sh');

if terminal_output ~= 0
    error('Error when compiling ARMOUR real time planner code!');
else
    fprintf("Successfully compiled ARMOUR\n\n");
end

cd ../../

%% Initialize real-time armtd planner for comparison
fprintf("Start compiling ARMTD\n\n");

cd kinova_simulator_interfaces/kinova_planner_realtime_armtd_comparison

armtd_buffer_path = [pwd, '/BufferPath.h'];
fid = fopen(armtd_buffer_path, 'w');
fprintf(fid, '#include <string>\n');
fprintf(fid, 'const std::string pathname = "%s/buffer/";', pwd);
fclose(fid);

terminal_output = system('./compile.sh');

if terminal_output ~= 0
    error('Error when compiling ARMTD real time planner code!');
else
    fprintf("Successfully compiled ARMTD\n\n");
end

cd ../../

%% Save current folder path as .mat file for other scripts to read
kinova_test_folder_path = pwd;

save('kinova_test_folder_path.mat', 'kinova_test_folder_path');

%% Create dir
mkdir kinova_simulator_interfaces/kinova_planner_realtime/buffer/
addpath kinova_simulator_interfaces/kinova_planner_realtime/buffer/
mkdir kinova_simulator_interfaces/kinova_planner_realtime_armtd_comparison/buffer/
addpath kinova_simulator_interfaces/kinova_planner_realtime_armtd_comparison/buffer/
