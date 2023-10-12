%% user parameters
close all; clear; clc;

file_location = '../results/hard_SL/' ;
scene_idx = 2;

verbosity = 0 ;
dimension = 3 ;

plot_start_and_end_config_only = true; % otherwise, animate trial.

%% automated from here
summary_files = dir([file_location, '/trial_*']);
load(summary_files(scene_idx).name);
    
%% plotting
figure(1) ; clf ; axis equal ; xlim([-1 1]); ylim([-1 1]); zlim([0 1.5]); grid on; hold on ; view(3) ;

A.animation_gif_filename = [file_location,'/','animation_',num2str(scene_idx),'.gif'];
frame_rate = A.animation_time_discretization / A.animation_playback_rate ;
start_gif = true;
W.plot();
A.plot_at_time(A.time(1));
breakpoint = 1;
for i = 1:10:(length(A.time)-1)
    A.plot_at_time(A.time(i));
   
    r_id = floor(i / 50) + 1;
%     for j = 1:7
%         hd{j} = trisurf(convhulln(P.info.sliced_FO_zono{r_id}{j}), ...
%                         P.info.sliced_FO_zono{r_id}{j}(:,1), ...
%                         P.info.sliced_FO_zono{r_id}{j}(:,2), ...
%                         P.info.sliced_FO_zono{r_id}{j}(:,3), ...
%                         'FaceColor',[0,0.5,0.5], ...
%                         'FaceAlpha',0.1, ...
%                         'EdgeColor',[0,0.5,0.5], ...
%                         'EdgeAlpha',0.2);
%     end

    fh = get(groot,'CurrentFigure') ;
    frame = getframe(fh) ;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if start_gif
        imwrite(imind,cm,A.animation_gif_filename,'gif', 'Loopcount',inf,...
                'DelayTime',frame_rate) ;
        start_gif = false ;
    else 
        imwrite(imind,cm,A.animation_gif_filename,'gif','WriteMode','append',...
                'DelayTime',frame_rate) ;
    end
    
%     for j = 1:7
%         delete(hd{j});
%     end
end