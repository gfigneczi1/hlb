function [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareData(segment_m, indexes, PARAMS)
     % input array
     input = [segment_m(:, indexes.OncomingTrafficType), ...
                segment_m(:, indexes.FrontTrafficType), ...
                segment_m(:, indexes.VelocityX), ...
                movmean(segment_m(:, indexes.AccelerationX),20), ...
                movmean(segment_m(:, indexes.YawRate),20), ...
                movmean(segment_m(:, indexes.LaneCurvature), 20), ...
                movmean(segment_m(:, indexes.c3), 200)];

         input(:,1) = movmean(input(:,1),50);
         input(:,2) = movmean(input(:,2),50);
    
    % output array
    output = zeros(size(input,1),numel(PARAMS.OUTPUT_SHIFT));
    delta = -segment_m(:,indexes.c0);
    for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)
        dT = mean(diff(segment_m(:, indexes.Relative_time)));
        dx = segment_m(:, indexes.VelocityX)*dT;        
        shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.OUTPUT_SHIFT(shiftID)./dx);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                output(shiftIDonOutput,shiftID) = delta(shiftOnOutput(shiftIDonOutput));
            end
        end
    end
    
    inputRaw = input;
    outputRaw = output;
    
    % SHUFFLE
    N = size(input,1);
    shuffledIndeces = randperm(N);
    input = input(shuffledIndeces,:);
    output = output(shuffledIndeces, :);

    % LIMIT DATA IF NEEDED
    % this is done before norm and central
    input = input(1:min(size(input,1),PARAMS.MAXIMUM_INPUT_SIZE),:);
    output = output(1:min(size(input,1),PARAMS.MAXIMUM_INPUT_SIZE), :);

    % NORM AND CENTRAL
    [output, c_out, s_out] = normAndCentral(output);
    [input, c_in, s_in] = normAndCentral(input);
end

function [input,c,s] = normAndCentral(input)
    for i=1:size(input,2)
         [input(:,i), c(i), s(i)] = normalize(input(:,i));
    end
end

