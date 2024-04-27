function [localization, road_type] = localization(rawData, tempFolder)
%LOCALIZATION Summary of this function goes here
%   Detailed explanation goes here
% root: path of the temp folder
root = tempFolder(1:end-6);
path = fullfile(root,'_matlab_evaluation','HLB','Inputs','reference_polygon.mat');
reference_polygon = load(path);
X_abs = rawData.LongPos_abs' * 40075000 .* cos(rawData.LatPos_abs'*pi()/180) / 360;
Y_abs = rawData.LatPos_abs' * 111.32*1000;
road_names = fieldnames(reference_polygon);
in_polygon_all = false(length(X_abs), 1);
road_type = zeros(1, length(X_abs));
for i=1:numel(road_names)
    % labeling roads, 1 for country road, 2 for highway
    if contains(road_names{i}, "r")
        road_type_name = 1;
    elseif contains(road_names{i}, "M")
        road_type_name = 2;
    elseif contains(road_names{i}, "ZZ")
        road_type_name = 3;
    elseif contains(road_names{i}, "ref")
        road_type_name = 4;
    end
    road_reference_polygon = reference_polygon.(road_names{i});
    in_polygon = inpolygon(X_abs, Y_abs, road_reference_polygon.X_poly,...
        road_reference_polygon.Y_poly);
    road_ind = find(in_polygon);
    if ~isempty(road_ind)
        road_type(road_ind) = road_type_name;
    end
    in_polygon_all = in_polygon_all | in_polygon;
end
localization = in_polygon_all';
end

