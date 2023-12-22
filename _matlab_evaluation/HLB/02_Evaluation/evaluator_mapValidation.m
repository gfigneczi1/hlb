function evaluator_mapValidation (segments,config)
    % Checking input data validity
    if (length(fieldnames(segments)) == 2 && segments.segments(1).signalStatus(2) == 1)
        transformed_coefficients = lane_model_validation (segments.validator, segments.segments(1).segment);

        segments.validated = segments.segments(1).segment;
        segments.validated.c01 = transformed_coefficients(:,1);
        segments.validated.c1 = transformed_coefficients(:,2);
        segments.validated.c2 = transformed_coefficients(:,3);
        segments.segments(1).segment.c01 = 0.5*(segments.segments(1).segment.c01_left+segments.segments(1).segment.c01_right);
        % Plotting the coefficient validations
        evaluate_plot_coefficient_validation(segments, config);
        % Calling the driver model
        [GT_U_ref, GT_Y_ref] = driverModelAndPlanner(config,segments.segment, "reference");
        [GT_U_val, GT_Y_val] = driverModelAndPlanner(config,segments.validated, "backValidation");

        % calculating KPIs
        calculate_KPIs(segments,config, GT_U_ref, GT_U_val, GT_Y_ref, GT_Y_val);
    elseif (isfield(segments, 'validator') && isfield(segments, 'segments'))
        % more files in the segments, one validator: producing missing
        % signals
        for i=1:length(segments.segments)
            if (segments.segments(i).signalStatus(1) == 1)
                transformed_coefficients = lane_model_validation_fineGrid (segments.validator, segments.segments(i).segment);
                segments.segments(i).segment.c01_left = transformed_coefficients(:,1)+1.875;
                segments.segments(i).segment.c01_right = transformed_coefficients(:,1)-1.875;
                segments.segments(i).segment.c1 = transformed_coefficients(:,2);
                segments.segments(i).segment.c2 = transformed_coefficients(:,3);
                segments.segments(i).segment.SteeringTorque = zeros(size(transformed_coefficients,1),1);
                segments.segments(i).segment.Left_Index = zeros(size(transformed_coefficients,1),1);
                segments.segments(i).segment.Right_Index = zeros(size(transformed_coefficients,1),1);
                % saving mat file
                segment = segments.segments(i).segment;
                segment = table2struct(segment,"ToScalar",true);
                signals = fieldnames(segment);
                for j=1:length(signals)
                    fieldname = convertCharsToStrings(signals(j));
                    segment.(fieldname) = segment.(fieldname)';
                end
                save (fullfile(config.root,strcat(segments.segments(i).name,"_extended.mat")),'segment');
            end
        end
    else
        disp('No validation possible, as the measurement file numbers are not valid');
    end
end

function calculate_KPIs (segments, config, GT_U_ref, GT_U_val, GT_Y_ref, GT_Y_val)
    % Calculating KPI-s
    results.GT_Y.avg = [mean(abs(GT_Y_ref(1,:)-GT_Y_val(1,:))); ...
        mean(abs(GT_Y_ref(2,:)-GT_Y_val(2,:)));
        mean(abs(GT_Y_ref(3,:)-GT_Y_val(3,:)));
        mean(abs(GT_Y_ref(4,:)-GT_Y_val(4,:)))];
    results.GT_Y.std = [std(abs(GT_Y_ref(1,:)-GT_Y_val(1,:)));
        std(abs(GT_Y_ref(2,:)-GT_Y_val(2,:)));
        std(abs(GT_Y_ref(3,:)-GT_Y_val(3,:)));
        std(abs(GT_Y_ref(4,:)-GT_Y_val(4,:)))];
    results.GT_U.avg = [mean(abs(GT_U_ref(1,:)-GT_U_val(1,:)))*1e-3; ...
        mean(abs(GT_U_ref(2,:)-GT_U_val(2,:)))*1e-3;
        mean(abs(GT_U_ref(3,:)-GT_U_val(3,:)))*1e-3];
    results.GT_U.std = [std(abs(GT_U_ref(1,:)-GT_U_val(1,:)))*1e-3;
        std(abs(GT_U_ref(2,:)-GT_U_val(2,:)))*1e-3;
        std(abs(GT_U_ref(3,:)-GT_U_val(3,:)))*1e-3];
    results.c0.avg = mean(abs(segments.validated.c01 - segments.segment.c01));
    results.c0.std = std(abs(segments.validated.c01 - segments.segment.c01));
    results.c2.avg = mean(abs(segments.validated.c2 - segments.segment.c2));
    results.c2.std = std(abs(segments.validated.c2 - segments.segment.c2));
    
    %7x3 curvature gradient model
    curveIndeces = abs(((GT_U_val(1,:)+GT_U_val(2,:)+GT_U_val(3,:))/3)) > (2.5e-4/1e-3);
    p1 = inv(GT_U_val(:,curveIndeces)*GT_U_val(:,curveIndeces)')*GT_U_val(:,curveIndeces)*(GT_Y_val(2,curveIndeces))';
    p2 = inv(GT_U_val(:,curveIndeces)*GT_U_val(:,curveIndeces)')*GT_U_val(:,curveIndeces)*(GT_Y_val(3,curveIndeces))';
    p3 = inv(GT_U_val(:,curveIndeces)*GT_U_val(:,curveIndeces)')*GT_U_val(:,curveIndeces)*(GT_Y_val(4,curveIndeces))';
    results.Poptim_val = [p1'; p2'; p3']; 
    
    curveIndeces = abs(((GT_U_ref(1,:)+GT_U_ref(2,:)+GT_U_ref(3,:))/3)) > (2.5e-4/1e-3);
    p1 = inv(GT_U_ref(:,curveIndeces)*GT_U_ref(:,curveIndeces)')*GT_U_ref(:,curveIndeces)*(GT_Y_ref(2,curveIndeces))';
    p2 = inv(GT_U_ref(:,curveIndeces)*GT_U_ref(:,curveIndeces)')*GT_U_ref(:,curveIndeces)*(GT_Y_ref(3,curveIndeces))';
    p3 = inv(GT_U_ref(:,curveIndeces)*GT_U_ref(:,curveIndeces)')*GT_U_ref(:,curveIndeces)*(GT_Y_ref(4,curveIndeces))';
    results.Poptim_ref = [p1'; p2'; p3']; 
    
    % plotting the surface
    plotPOptim_service(results.Poptim_ref, results.Poptim_val, config);
    
    disp('KPI-s for map validation:');
    disp('GT_Y - avg:');
    disp(results.GT_Y.avg);
    disp('GT_Y - std:');
    disp(results.GT_Y.std);
    disp('GT_U kappa - avg:');
    disp(results.GT_U.avg);
    disp('GT_U kappa - std:');
    disp(results.GT_U.std);
    disp('c0 - avg:');
    disp(results.c0.avg);
    disp('c0 - std:');
    disp(results.c0.std);
    disp('c2 - avg:');
    disp(results.c2.avg);
    disp('c2 - std:');
    disp(results.c2.std);
    disp('Poptim - ref:');
    disp(results.Poptim_ref);
    disp('Poptim - val:');
    disp(results.Poptim_val);
    
    save(fullfile(config.root,'plots','map_validation_results.mat'),'results');
end

function [transformed_coefficients] = lane_model_validation (validator, segment)
%Inputs:
% - reference: array of N x 8, namely: t X Y theta c0 c1 c2 x0
% - validator: array of N x 8, namely: t X Y theta [0]-1x3 [0]
% where N is the number of samples


X_abs_val = validator.segment.LongPos_abs * 40075000 .* cos(validator.segment.LatPos_abs*pi()/180) / 360;
Y_abs_val = validator.segment.LatPos_abs * 111.32*1000;
theta_calc_val = atan(tan(validator.segment.theta_calc));

X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
Y_abs = segment.LatPos_abs * 111.32*1000;
theta_calc = atan(tan(segment.theta_calc));

prepareInputForPlanner(validator.segment);

for i=1:size(segment,1)
    client = [X_abs(i) Y_abs(i) theta_calc(i)];
    % finding closest point
    distances = ((X_abs_val-client(1)).^2+(Y_abs_val-client(2)).^2).^0.5;
    if (min(distances) < 10)
        closest_point_idx = find(distances==min(distances),1);
        base = [X_abs_val(closest_point_idx(1,1)) Y_abs_val(closest_point_idx(1,1)) theta_calc_val(closest_point_idx(1,1))];
        dtheta_in = 0; % change to pi() when the vehicles are in opposite direction
        variation = "validation";
        % change to c02_left from c01_right when vehicles are in opposite
        % direction
        coefficients = [0.5*(validator.segment.c01_left(closest_point_idx(1,1))+validator.segment.c01_right(closest_point_idx(1,1))) validator.segment.c1(closest_point_idx(1,1)) validator.segment.c2(closest_point_idx(1,1))];
        transformed_coefficients(i,1:3) = lane_border_interpolation_2ndorder(base,client,coefficients, dtheta_in, variation);
    else
        transformed_coefficients(i,1:3) = [nan nan nan];
    end
end


        
end

function [transformed_coefficients] = lane_model_validation_fineGrid (validator, segment)
% using splin to create a fine grid on the coefficients from the validator.
% then, we search for the closest point from the segment, and that is
% to calculate the coefficients, the following way:
% c0: closest point is searched, then transformation based on validator ego
% pose is executed.
% c1: only validator theta is considered
% c2: no transformation needed.
validator.segment.X_abs = validator.segment.LongPos_abs * 40075000 .* cos(validator.segment.LatPos_abs*pi()/180) / 360;
validator.segment.Y_abs = validator.segment.LatPos_abs * 111.32*1000;

X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
Y_abs = segment.LatPos_abs * 111.32*1000;
segment.X_abs = X_abs;
segment.Y_abs = Y_abs;

% spline fitting is only possibly if X_abs coordinates are monotonously
% increasing
% checking direction:
dir = sign(mean(diff(validator.segment.X_abs)));
if (dir == 1)
    validatorNotMon = numel(find(diff(validator.segment.X_abs)<=0)) > 0;
else
    validatorNotMon = numel(find(diff(validator.segment.X_abs)>=0)) > 0;
end
dir = sign(mean(diff(segment.X_abs)));
if (dir == 1)
    segmentNotMon = numel(find(diff(segment.X_abs)<=0)) > 0;
else
    segmentNotMon = numel(find(diff(segment.X_abs)>=0)) > 0;
end
if (~validatorNotMon && ~segmentNotMon)
    disp('Using the fine grid based new solution for coefficient transformation');
    yy = spline(validator.segment.X_abs, validator.segment.Y_abs, segment.X_abs);
    tt = spline(validator.segment.X_abs, validator.segment.theta_calc,segment.X_abs);
    c0_validator = 0.5*(validator.segment.c01_left+validator.segment.c01_right);
    c1_validator = validator.segment.c1;
    c2_validator = validator.segment.c2;
    c00 = spline(validator.segment.X_abs,c0_validator,segment.X_abs);
    c11 = spline(validator.segment.X_abs,c1_validator,segment.X_abs);
    c22 = spline(validator.segment.X_abs,c2_validator,segment.X_abs);
    dy = (segment.Y_abs - yy).*cos(tt);
    dtheta = segment.theta_calc - tt;
    transformed_coefficients(:,1) = c00- dy./cos(dtheta);
    transformed_coefficients(:,2) = c11 - dtheta;
    transformed_coefficients(:,3) = c22;          
else
    disp('Coefficient transformation with spline is not possible, as X values are not monotonic');
    disp('Using the old solution instead');
    transformed_coefficients = lane_model_validation (validator, segment);
end


end

function plotPOptim_service (P_ref, P_val, config)

    f = figure();
    title('P surface','FontSize',18);
    subplot(1,2,1);
    surf(P_ref);
    title('P_ref');
    subplot(1,2,2);
    surf(P_val);
    title('P_val');

    savefig(f, fullfile(config.root, 'plots','P_surf.fig'));
    saveas(f, fullfile(config.root, 'plots','P_surf.png'));

    close(f);

end