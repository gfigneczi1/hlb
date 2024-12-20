function [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareDataPHTPM(segment_m, indexes, PARAMS)
     % input array
     ladVector = PARAMS.Lp*segment_m(:,indexes.VelocityX);
            ye = segment_m(:,indexes.c0)+ ...
                tan(segment_m(:,indexes.LaneOrientation)).*ladVector + ...
                segment_m(:,indexes.LaneCurvature)/2.*ladVector.^2 + ...
                segment_m(:,indexes.c3)/6.*ladVector.^3;
            thetaFP = atan(ye/PARAMS.Lp);

            determinant = 4*(segment_m(:,indexes.LaneCurvature)/2).^2 - 12*tan(segment_m(:,indexes.LaneOrientation)).*segment_m(:,indexes.c3)/6;
            determinant(determinant<0) = nan;
    
            x1 = (-2*segment_m(:,indexes.LaneCurvature)/2 + determinant)./(segment_m(:,indexes.c3));
            x2 = (-2*segment_m(:,indexes.LaneCurvature)/2 - determinant)./(segment_m(:,indexes.c3));
            x = min(max(0,x1),max(0,x2));
            x(isnan(determinant)) = PARAMS.Lp*segment_m(isnan(determinant), indexes.VelocityX);
            x(x==0) = PARAMS.Lp*segment_m(x==0, indexes.VelocityX);
    
            kappa_r = segment_m(:,indexes.LaneCurvature)+segment_m(:,indexes.c3).*x;
            thetaTP = atan(x.*kappa_r);

            input = [movmean(thetaTP,100), ...
                    movmean(thetaFP,100), ...
                    segment_m(:, indexes.VelocityX), ...
                    movmean(segment_m(:, indexes.Acceleration_Y),100), ...
                    movmean(segment_m(:, indexes.YawRate),100), ...
                    movmean(-segment_m(:, indexes.LaneOrientation), 20)];
    
    % output array
    output = -0.5*(segment_m(:,indexes.LaneEdgePositionLeft)+segment_m(:,indexes.LaneEdgePositionRight));
    
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

