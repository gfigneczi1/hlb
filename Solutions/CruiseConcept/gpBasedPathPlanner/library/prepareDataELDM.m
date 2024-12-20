function [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareDataELDM(segment_m, indexes, PARAMS)
     % input array - based on PARAMS.npDistances
     dT = mean(diff(segment_m(:, indexes.Relative_time)));
     dx = segment_m(:, indexes.VelocityX)*dT;   
     filteredCurvature = movmean(segment_m(:, indexes.LaneCurvature), 20);
     for np=length(PARAMS.npDistances):-1:1
         % calculate the average curvature between the node points
         if (np==1)
             shiftOnOutput_1 = [1:1:size(segment_m(:,indexes.Relative_time),1)]';
             shiftOnOutput_2 = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.npDistances(np)./dx);
         else
             shiftOnOutput_1 = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.npDistances(np-1)./dx);
             shiftOnOutput_2 = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.npDistances(np)./dx);
         end
         for i=1:length(filteredCurvature)
             if (shiftOnOutput_2(i) > length(filteredCurvature))
                 break;
             end
             input(i,np) = mean(filteredCurvature(shiftOnOutput_1(i):shiftOnOutput_2(i)));
         end
         % last n points missing due to missing preview information 
     end
  
    % output array
    output = zeros(size(input,1),numel(PARAMS.npDistances));
    delta = -segment_m(:,indexes.c0);
    for shiftID=1:numel(PARAMS.npDistances)
           
        shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.npDistances(shiftID)./dx);
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

