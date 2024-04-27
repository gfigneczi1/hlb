clear; close all;

load('C:\database\KDP_HLB_GP\Driver_with_traffic.mat');
config.root = "./";

for i=1:length(segments.segments)
    segment = segments.segments(i).segment;
    name = segments.segments(i).name;

    % transforming the struct to array for lower calculation time
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    
    temp_folder_path = config.root;
    plots_folder_name = 'plots';
    set(0,'DefaultFigureVisible','off');

    [traj, segments] = trajectory_planner(anchorPoints, indeces, ref, ver)

end

