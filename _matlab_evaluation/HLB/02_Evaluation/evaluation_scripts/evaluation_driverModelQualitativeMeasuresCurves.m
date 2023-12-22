function [curves] = evaluation_driverModelQualitativeMeasures(curve,  name_fig, name_png)

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-07

displacement = curve.displacement;
T = [cos(curve.orient(1)) sin(curve.orient(1)); -sin(curve.orient(1)) cos(curve.orient(1))];
ref = [curve.ref(:,1)-curve.cor(1,1) curve.ref(:,2)-curve.cor(1,2)]*T';
cor = [curve.cor(:,1)-curve.cor(1,1) curve.cor(:,2)-curve.cor(1,2)]*T';
traj = [curve.traj(:,1)-curve.cor(1,1) curve.traj(:,2)-curve.cor(1,2)]*T';
orient = curve.orient - curve.orient(1);

% KPI-1: side correctness
refSignedDeviation = sign(ref(:,2)-cor(:,2)).*(((ref(:,1)-cor(:,1)).^2+(ref(:,2)-cor(:,2)).^2).^0.5);
trajSignedDeviation = sign(traj(:,2)-cor(:,2)).*(((traj(:,1)-cor(:,1)).^2+(traj(:,2)-cor(:,2)).^2).^0.5);

temp_folder_path = fullfile('..','..','_temp');
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

T = [cos(curve.orient(1)) sin(curve.orient(1)); -sin(curve.orient(1)) cos(curve.orient(1))];
corLocal = (curve.cor-curve.cor(1,:))*T';
corLeftLocal = (curve.corLeft-curve.cor(1,:))*T';
corRightLocal = (curve.corRight-curve.cor(1,:))*T';

% Calculating side correctness
deadzone = (abs(refSignedDeviation) < 0.05); %+(abs(curve.curv)<2e-5);
deadzoneTraj = (abs(trajSignedDeviation) < 0.05); %+(abs(curve.curv)<2e-5);
% term of actual side correctness outside of deadzone
correctnessVector = sign(refSignedDeviation(deadzone==0)).*sign(trajSignedDeviation(deadzone==0));
% term of accuracy inside the deadzone
accuracyVector = deadzoneTraj(deadzone) > 0;
curves.sideCorrectness = (numel(find(correctnessVector==1)) + numel(find(accuracyVector==1)))/length(refSignedDeviation);

% KPI - 2: generic deviations
positionDeviation = ((curve.ref(:,1)-curve.traj(:,1)).^2 + (curve.ref(:,2)-curve.traj(:,2)).^2).^0.5;
curves.avgPositionDeviation = mean(positionDeviation);
curves.medianPositionDeviation = median(positionDeviation);
curves.maxPositionDeviation = max(positionDeviation);

% border violation
curves.borderViolation = numel(find((abs(trajSignedDeviation)>(1.875-0.85)).*(abs(abs(refSignedDeviation)) < 1.25))) / size(trajSignedDeviation,1);
if (curves.borderViolation == 0)
    curves.borderDistance = (1.875-0.85) - max(abs(trajSignedDeviation));
else
    curves.borderDistance = 0;
end
% resulting steering angle
refCurvature = -1*movmean([0; movmean([0; diff(diff(curve.ref(:,2))./diff(curve.ref(:,1)))]./diff(curve.ref(:,1)),250)], 50);
trajCurvature = -1*movmean([0; movmean([0; diff(diff(curve.traj(:,2))./diff(curve.traj(:,1)))]./diff(curve.traj(:,1)),250)], 50);
corCurvature = -1*movmean([0; movmean([0; diff(diff(curve.cor(:,2))./diff(curve.cor(:,1)))]./diff(curve.cor(:,1)),250)],50);
curves.nominalCurvature = mean(abs(corCurvature));
curves.avgCurvatureDeviation = mean(abs(refCurvature-trajCurvature));
curves.integratedCurvatureDeviation = trapz((refCurvature-trajCurvature).^2)/length(refCurvature);
curves.avgCurvatureDeviationRefToCor = mean(abs(refCurvature-corCurvature));
curves.curveLength = sum((diff(cor(:,1)).^2+diff(cor(:,2)).^2).^0.5);

set(0,'DefaultFigureVisible','off');

end

