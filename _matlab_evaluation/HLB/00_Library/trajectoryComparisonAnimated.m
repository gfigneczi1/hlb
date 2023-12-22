function trajectoryComparisonAnimated(LDM_out,LDMPoly_out,TRC_in, model_out)
%TRAJECTORYCOMPARISONANIMATED Summary of this function goes here
%   Detailed explanation goes here
% animation time
startIdx = 1300;
endIdx = 5500;

TRC_PLOT = false;

lengths = [size(LDM_out.trajectoryEgoFrame.Data,3), size(TRC_in.length,2)];
endIdxMax = min(lengths);

LDM_out.trajectoryEgoFrame.Data = LDM_out.trajectoryEgoFrame.Data(:,:,1:endIdxMax);
LDMPoly_out.trajectoryEgoFrame.Data = LDMPoly_out.trajectoryEgoFrame.Data(:,:,1:endIdxMax);
LDMPoly_out.trajectoryCoefficients.Data = LDMPoly_out.trajectoryCoefficients.Data(1:endIdxMax,:);
model_out.egoPoseGlobalFrame.Data = model_out.egoPoseGlobalFrame.Data(:,:,1:endIdxMax);


% downsampling from 1 ms to 10 ms
TRC_in.FirstSegment_c0 = TRC_in.FirstSegment_c0(1:10:endIdxMax*10);
TRC_in.FirstSegment_c1 = TRC_in.FirstSegment_c1(1:10:endIdxMax*10);
TRC_in.FirstSegment_c2 = TRC_in.FirstSegment_c2(1:10:endIdxMax*10);
TRC_in.FirstSegment_c3 = TRC_in.FirstSegment_c3(1:10:endIdxMax*10);
TRC_in.FirstSegment_c4 = TRC_in.FirstSegment_c4(1:10:endIdxMax*10);
TRC_in.FirstSegment_c5 = TRC_in.FirstSegment_c5(1:10:endIdxMax*10);

TRC_in.SecondSegment_c0 = TRC_in.SecondSegment_c0(1:10:endIdxMax*10);
TRC_in.SecondSegment_c1 = TRC_in.SecondSegment_c1(1:10:endIdxMax*10);
TRC_in.SecondSegment_c2 = TRC_in.SecondSegment_c2(1:10:endIdxMax*10);
TRC_in.SecondSegment_c3 = TRC_in.SecondSegment_c3(1:10:endIdxMax*10);

TRC_in.breakPoint = TRC_in.breakPoint(1:10:endIdxMax*10);
TRC_in.length = TRC_in.length(1:10:endIdxMax*10);
TRC_in.q_T0 = TRC_in.q_T0(1:10:endIdxMax*10);

%TRC_traj = LDM_out.trajectoryEgoFrame.Data;
LDM_fitted_poly = LDM_out.trajectoryEgoFrame.Data;

%for i=1:size(TRC_traj(:,1,:),3)
%    TRC_traj(1:TRC_traj(end,1,i),2,i) = TRC_in.SecondSegment_c0(i)' + TRC_traj(1:TRC_traj(end,1,i),1,i) .* TRC_in.SecondSegment_c1(i)' + (TRC_traj(1:TRC_traj(end,1,i),1,i).^2) .* (TRC_in.SecondSegment_c2(i)' .* 1/2) + (TRC_traj(1:TRC_traj(end,1,i),1,i).^3) .* (TRC_in.SecondSegment_c3(i)' .* 1/6);
%end
for i=1:size(LDM_fitted_poly(:,1,:),3)
    LDM_fitted_poly(1:LDM_fitted_poly(end,1,i),2,i) = LDM_out.trajectoryFittedCoefficients.Data(1,1,i)' + LDM_fitted_poly(1:LDM_fitted_poly(end,1,i),1,i) .* LDM_out.trajectoryFittedCoefficients.Data(1,2,i)' + (LDM_fitted_poly(1:LDM_fitted_poly(end,1,i),1,i).^2) .* (LDM_out.trajectoryFittedCoefficients.Data(1,3,i)') + (LDM_fitted_poly(1:LDM_fitted_poly(end,1,i),1,i).^3) .* (LDM_out.trajectoryFittedCoefficients.Data(1,4,i)');
end

figure(1); hold on; grid on;
for n=startIdx:2:endIdx
    plot((LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),1,n).*cos(model_out.egoPoseGlobalFrame.Data(3,1,n)) - sin(model_out.egoPoseGlobalFrame.Data(3,1,n)).*LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(1,1,n), (LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),1,n).*sin(model_out.egoPoseGlobalFrame.Data(3,1,n)) + cos(model_out.egoPoseGlobalFrame.Data(3,1,n)).*LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(2,1,n), "bo", "LineWidth", 2);
    hold on;
    %plot((LDMPoly_out.trajectoryEgoFrame.Data(1:LDMPoly_out.trajectoryEgoFrame.Data(301,1,n),1,n).*cos(model_out.egoPoseGlobalFrame.Data(3,1,n)) - sin(model_out.egoPoseGlobalFrame.Data(3,1,n)).*LDMPoly_out.trajectoryEgoFrame.Data(1:LDMPoly_out.trajectoryEgoFrame.Data(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(1,1,n), (LDMPoly_out.trajectoryEgoFrame.Data(1:LDMPoly_out.trajectoryEgoFrame.Data(301,1,n),1,n).*sin(model_out.egoPoseGlobalFrame.Data(3,1,n)) + cos(model_out.egoPoseGlobalFrame.Data(3,1,n)).*LDMPoly_out.trajectoryEgoFrame.Data(1:LDMPoly_out.trajectoryEgoFrame.Data(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(2,1,n), "ro", "LineWidth", 1);
    %hold on;
    plot((LDM_fitted_poly(1:LDM_fitted_poly(301,1,n),1,n).*cos(model_out.egoPoseGlobalFrame.Data(3,1,n)) - sin(model_out.egoPoseGlobalFrame.Data(3,1,n)).*LDM_fitted_poly(1:LDM_fitted_poly(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(1,1,n), (LDM_fitted_poly(1:LDM_fitted_poly(301,1,n),1,n).*sin(model_out.egoPoseGlobalFrame.Data(3,1,n)) + cos(model_out.egoPoseGlobalFrame.Data(3,1,n)).*LDM_fitted_poly(1:LDM_fitted_poly(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(2,1,n), "co", "LineWidth", 1);
    hold on;
    if TRC_PLOT == true
        plot((TRC_traj(1:TRC_traj(301,1,n),1,n).*cos(model_out.egoPoseGlobalFrame.Data(3,1,n)) - sin(model_out.egoPoseGlobalFrame.Data(3,1,n)).*TRC_traj(1:TRC_traj(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(1,1,n), (TRC_traj(1:TRC_traj(301,1,n),1,n).*sin(model_out.egoPoseGlobalFrame.Data(3,1,n)) + cos(model_out.egoPoseGlobalFrame.Data(3,1,n)).*TRC_traj(1:TRC_traj(301,1,n),2,n)) + model_out.egoPoseGlobalFrame.Data(2,1,n), "cx-", "LineWidth", 1);
        hold on;
        plot((TRC_traj(1:27,1,n).*cos(model_out.egoPoseGlobalFrame.Data(3,1,n)) - sin(model_out.egoPoseGlobalFrame.Data(3,1,n)).*TRC_traj(1:27,2,n)) + model_out.egoPoseGlobalFrame.Data(1,1,n), (TRC_traj(1:27,1,n).*sin(model_out.egoPoseGlobalFrame.Data(3,1,n)) + cos(model_out.egoPoseGlobalFrame.Data(3,1,n)).*TRC_traj(1:27,2,n)) + model_out.egoPoseGlobalFrame.Data(2,1,n), "kx-", "LineWidth", 2);
        hold on;
    end
    plot(model_out.egoPoseGlobalFrame.Data(1,1,n),model_out.egoPoseGlobalFrame.Data(2,1,n),"gx", "MarkerSize", 15, "LineWidth", 4); hold on;
    if TRC_PLOT == true
        legend('clothoid ', '3rdpoly', 'TRCpoly', 'TRC 25m','ego');
    else
        legend('clothoid ', 'fitted poly','ego');
    end
    %distances_mean = mean( ( (LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),1,n) - LDMPoly_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),1,n)).^2 + (LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,2,n),2,n) - LDMPoly_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,2,n),2,n)).^2 ).^0.5 );
    %distances_max = max( ( (LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),1,n) - LDMPoly_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,1,n),1,n)).^2 + (LDM_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,2,n),2,n) - LDMPoly_out.trajectoryEgoFrame.Data(1:LDM_out.trajectoryEgoFrame.Data(301,2,n),2,n)).^2 ).^0.5 );
    %annotation('textbox',[.35 .7 .2 .2],'String',strcat('mean error: ', num2str(round(distances_mean*100,2)), ' cm'),'FitBoxToText','on');
    %annotation('textbox',[.35 .66 .2 .2],'String',strcat('max error: ', num2str(round(distances_max*100,2)), ' cm'),'FitBoxToText','on');
    grid on;
    pause(0.02);
    clf;
    disp(n);
end

end

