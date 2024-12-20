function [segment_m, indexes, expectedOutput] = conditionData(modelID, segment_m, indexes, PARAMS)

% newly added channels are set to default
indexes.LDM_N = 0;
indexes.LDM_1 = 0;

switch modelID
    case "sgp"
        [~, ~, inputRaw, expectedOutput] = prepareData(segment_m, indexes, PARAMS);
        segment_m(:,indexes.AccelerationX) = inputRaw(:,4);
        segment_m(:,indexes.VelocityX) = inputRaw(:,3);
        segment_m(:,indexes.YawRate) = inputRaw(:,5);
        segment_m(:,indexes.OncomingTrafficType) = inputRaw(:,1);
        segment_m(:,indexes.FrontTrafficType) = inputRaw(:,2);
        segment_m(:,indexes.LaneCurvature) = inputRaw(:,6);
        segment_m(:,indexes.c3) = inputRaw(:,7);
    case "narx"
        [~, ~, inputRaw, expectedOutput] = prepareDataNARMAX(segment_m, indexes, PARAMS);
        segment_m(:,indexes.AccelerationX) = inputRaw(:,4);
        segment_m(:,indexes.VelocityX) = inputRaw(:,3);
        segment_m(:,indexes.YawRate) = inputRaw(:,5);
        segment_m(:,indexes.OncomingTrafficType) = inputRaw(:,1);
        segment_m(:,indexes.FrontTrafficType) = inputRaw(:,2);
        segment_m(:,indexes.LaneCurvature) = inputRaw(:,6);
        segment_m(:,indexes.c3) = inputRaw(:,7);
    case "arx"
        [~, ~, inputRaw, expectedOutput] = prepareDataNARMAX(segment_m, indexes, PARAMS);
        segment_m(:,indexes.AccelerationX) = inputRaw(:,4);
        segment_m(:,indexes.VelocityX) = inputRaw(:,3);
        segment_m(:,indexes.YawRate) = inputRaw(:,5);
        segment_m(:,indexes.OncomingTrafficType) = inputRaw(:,1);
        segment_m(:,indexes.FrontTrafficType) = inputRaw(:,2);
        segment_m(:,indexes.LaneCurvature) = inputRaw(:,6);
        segment_m(:,indexes.c3) = inputRaw(:,7);
    case "phtpm"
        [~, ~, inputRaw, expectedOutput] = prepareDataPHTPM(segment_m, indexes, PARAMS);
        segment_m(:, indexes.VelocityX) = inputRaw(:,3);
        segment_m(:, indexes.Acceleration_Y) = inputRaw(:,4);
        segment_m(:,indexes.YawRate) = inputRaw(:,5);
        segment_m(:,indexes.LaneOrientation) = inputRaw(:,6);

        % thetaTP and thetaFP shall be calculated in real time, but for
        % reference the original data is added
        segment_m(:, end+1) = inputRaw(:,1);
        indexes.thetaTP = size(segment_m,2);
        segment_m(:, end+1) = inputRaw(:,2);
        indexes.thetaFP = size(segment_m,2);
    case "lrm"
        [~, ~, inputRaw, expectedOutput] = prepareDataNARMAX(segment_m, indexes, PARAMS);
        segment_m(:,indexes.AccelerationX) = inputRaw(:,4);
        segment_m(:,indexes.VelocityX) = inputRaw(:,3);
        segment_m(:,indexes.YawRate) = inputRaw(:,5);
        segment_m(:,indexes.OncomingTrafficType) = inputRaw(:,1);
        segment_m(:,indexes.FrontTrafficType) = inputRaw(:,2);
        segment_m(:,indexes.LaneCurvature) = inputRaw(:,6);
        segment_m(:,indexes.c3) = inputRaw(:,7);
    case "ldm"
        [~, ~, inputRaw, expectedOutput] = prepareDataELDM(segment_m, indexes, PARAMS);
        indexes.LDM_1 = size(segment_m,2)+1;
        segment_m = segment_m(1:size(inputRaw,1),:);
        for n=1:size(inputRaw,2)
            segment_m(:,end+1) = inputRaw(:,n);
        end
        indexes.LDM_N = size(segment_m,2);   
    case "eldm"
        [~, ~, inputRaw, expectedOutput] = prepareDataELDM(segment_m, indexes, PARAMS);
        indexes.LDM_1 = size(segment_m,2)+1;
        segment_m = segment_m(1:size(inputRaw,1),:);
        for n=1:size(inputRaw,2)
            segment_m(:,end+1) = inputRaw(:,n);
        end
        indexes.LDM_N = size(segment_m,2); 
end

end

