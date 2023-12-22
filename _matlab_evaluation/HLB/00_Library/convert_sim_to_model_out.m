function convert_sim_to_model_out(model_out, model_out_simulink)
%CONVERT_SIM_TO_MODEL_OUT Summary of this function goes here
%   Detailed explanation goes here
% load mat converted from D97
load('C:\Users\DOT1BP\Desktop\CM_sims_LDM_2023\LDM_P0\LDM_SFun_model_out.mat');
% planner traj
varDefaultPlannerX = 'trajectoryPlannerFrame_X_';
varDefaultPlannerY = 'trajectoryPlannerFrame_Y_';
% global traj
varDefaultGlobalX = 'trajectoryGlobalFrame_X_';
varDefaultGlobalY = 'trajectoryGlobalFrame_Y_';
% ego traj
varDefaultEgoX = 'trajectoryEgoFrame_X_';
varDefaultEgoY = 'trajectoryEgoFrame_Y_';
% planner corridor
varDefaultCPlannerX = 'corridorPlannerFrame_X_';
varDefaultCPlannerY = 'corridorPlannerFrame_Y_';
varDefaultCPlannerC1 = 'corridorPlannerFrame_c1_';
varDefaultCPlannerC2 = 'corridorPlannerFrame_c2_';
% global corridor
varDefaultCGlobalX = 'corridorGlobalFrame_X_';
varDefaultCGlobalY = 'corridorGlobalFrame_Y_';
varDefaultCGlobalC1 = 'corridorGlobalFrame_c1_';
varDefaultCGlobalC2 = 'corridorGlobalFrame_c2_';
% ego corridor
varDefaultCEgoX = 'corridorEgoFrame_X_';
varDefaultCEgoY = 'corridorEgoFrame_Y_';
varDefaultCEgoC1 = 'corridorEgoFrame_c1_';
varDefaultCEgoC2 = 'corridorEgoFrame_c2_';

model_out.trajectoryPlannerFrame.trajectoryPointsX = zeros(300,26301);
model_out.trajectoryPlannerFrame.trajectoryPointsY = zeros(300,26301);
model_out.trajectoryEgoFrame.trajectoryPointsX = zeros(300,26301);
model_out.trajectoryEgoFrame.trajectoryPointsY = zeros(300,26301);
model_out.trajectoryGlobalFrame.trajectoryPointsX = zeros(300,26301);
model_out.trajectoryGlobalFrame.trajectoryPointsY = zeros(300,26301);
model_out.corridorPlannerFrame.corridorPointsX = zeros(300,26301);
model_out.corridorPlannerFrame.corridorPointsY = zeros(300,26301);
model_out.corridorPlannerFrame.corridorC1 = zeros(300,26301);
model_out.corridorPlannerFrame.corridorC2 = zeros(300,26301);
model_out.corridorEgoFrame.corridorPointsX = zeros(300,26301);
model_out.corridorEgoFrame.corridorPointsY = zeros(300,26301);
model_out.corridorEgoFrame.corridorC1 = zeros(300,26301);
model_out.corridorEgoFrame.corridorC2 = zeros(300,26301);
model_out.corridorGlobalFrame.corridorPointsX = zeros(300,26301);
model_out.corridorGlobalFrame.corridorPointsY = zeros(300,26301);
model_out.corridorGlobalFrame.corridorC1 = zeros(300,26301);
model_out.corridorGlobalFrame.corridorC2 = zeros(300,26301);
for i=1:300
    %% trajectory
    varNamePlannerX = append(varDefaultPlannerX, num2str(i-1));
    varNamePlannerY = append(varDefaultPlannerY, num2str(i-1));
    varNameGlobalX = append(varDefaultGlobalX, num2str(i-1));
    varNameGlobalY = append(varDefaultGlobalY, num2str(i-1));
    varNameEgoX = append(varDefaultEgoX, num2str(i-1));
    varNameEgoY = append(varDefaultEgoY, num2str(i-1));
    eval(['model_out.trajectoryPlannerFrame.trajectoryPointsX(i,:) = ', varNamePlannerX, ';']);
    eval(['model_out.trajectoryPlannerFrame.trajectoryPointsY(i,:) = ', varNamePlannerY, ';']);
    eval(['model_out.trajectoryGlobalFrame.trajectoryPointsX(i,:) = ', varNameGlobalX, ';']);
    eval(['model_out.trajectoryGlobalFrame.trajectoryPointsY(i,:) = ', varNameGlobalY, ';']);
    eval(['model_out.trajectoryEgoFrame.trajectoryPointsX(i,:) = ', varNameEgoX, ';']);
    eval(['model_out.trajectoryEgoFrame.trajectoryPointsY(i,:) = ', varNameEgoY, ';']);
    %% corridor
    varNameCPlannerX = append(varDefaultCPlannerX, num2str(i-1));
    varNameCPlannerY = append(varDefaultCPlannerY, num2str(i-1));
    varNameCPlannerC1 = append(varDefaultCPlannerC1, num2str(i-1));
    varNameCPlannerC2 = append(varDefaultCPlannerC2, num2str(i-1));
    varNameCGlobalX = append(varDefaultCGlobalX, num2str(i-1));
    varNameCGlobalY = append(varDefaultCGlobalY, num2str(i-1));
    varNameCGlobalC1 = append(varDefaultCGlobalC1, num2str(i-1));
    varNameCGlobalC2 = append(varDefaultCGlobalC2, num2str(i-1));
    varNameCEgoX = append(varDefaultCEgoX, num2str(i-1));
    varNameCEgoY = append(varDefaultCEgoY, num2str(i-1));
    varNameCEgoC1 = append(varDefaultCEgoC1, num2str(i-1));
    varNameCEgoC2 = append(varDefaultCEgoC2, num2str(i-1));
    eval(['model_out.corridorPlannerFrame.corridorPointsX(i,:) = ', varNameCPlannerX, ';']);
    eval(['model_out.corridorPlannerFrame.corridorPointsY(i,:) = ', varNameCPlannerY, ';']);
    eval(['model_out.corridorPlannerFrame.corridorC1(i,:) = ', varNameCPlannerC1, ';']);
    eval(['model_out.corridorPlannerFrame.corridorC2(i,:) = ', varNameCPlannerC2, ';']);
    eval(['model_out.corridorGlobalFrame.corridorPointsX(i,:) = ', varNameCGlobalX, ';']);
    eval(['model_out.corridorGlobalFrame.corridorPointsY(i,:) = ', varNameCGlobalY, ';']);
    eval(['model_out.corridorGlobalFrame.corridorC1(i,:) = ', varNameCGlobalC1, ';']);
    eval(['model_out.corridorGlobalFrame.corridorC2(i,:) = ', varNameCGlobalC2, ';']);
    eval(['model_out.corridorEgoFrame.corridorPointsX(i,:) = ', varNameCEgoX, ';']);
    eval(['model_out.corridorEgoFrame.corridorPointsY(i,:) = ', varNameCEgoY, ';']);
    eval(['model_out.corridorEgoFrame.corridorC1(i,:) = ', varNameCEgoC1, ';']);
    eval(['model_out.corridorEgoFrame.corridorC2(i,:) = ', varNameCEgoC2, ';']);
%save('model_out.mat','-struct','model_out');
end
%% distance between the middle of the lane and trajectory
distances = mean(((model_out.corridorGlobalFrame.corridorPointsX(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4))),4) - model_out.trajectoryGlobalFrame.trajectoryPointsX(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4))),4)).^2 + (model_out.corridorGlobalFrame.corridorPointsY(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4))),4) - model_out.trajectoryGlobalFrame.trajectoryPointsY(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4))),4)).^2).^0.5);
for i=25:25:length(corridorEgoFrame_c1_0)
    distances = [distances, mean(((model_out.corridorGlobalFrame.corridorPointsX(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i))),i) - model_out.trajectoryGlobalFrame.trajectoryPointsX(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i))),i)).^2 + (model_out.corridorGlobalFrame.corridorPointsY(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i))),i) - model_out.trajectoryGlobalFrame.trajectoryPointsY(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i))),i)).^2).^0.5)];
end
mean_d = mean(distances(1:end));
max_d = max(distances(1:end));
min_d = min(distances(1:end));
std_d = std(distances(1:end));
disp(['Mean of the difference between reference lane and trajectory ', num2str(mean_d), ' meters']);
disp(['Maximum of the difference between reference lane and trajectory ', num2str(max_d), ' meters']);
disp(['Minimum of the difference between reference lane and trajectory ', num2str(min_d), ' meters']);
disp(['Standard deviation of the difference between reference lane and trajectory ', num2str(std_d), ' meters']);
%% distance between the S-function LDM traj vs. Simulink LDM traj
distances2 = mean(((model_out_simulink.trajectoryGlobalFrame.Data(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4)))-1,1,4) - model_out.trajectoryGlobalFrame.trajectoryPointsX(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4)))-1,4)).^2 + (model_out_simulink.trajectoryGlobalFrame.Data(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4)))-1,2,4) - model_out.trajectoryGlobalFrame.trajectoryPointsY(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,4)))-1,4)).^2).^0.5);
for i=25:25:length(corridorEgoFrame_c1_0)
    distances2 = [distances2, mean(((model_out_simulink.trajectoryGlobalFrame.Data(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i)))-1,1,i) - model_out.trajectoryGlobalFrame.trajectoryPointsX(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i)))-1,i)).^2 + (model_out_simulink.trajectoryGlobalFrame.Data(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i)))-1,2,i) - model_out.trajectoryGlobalFrame.trajectoryPointsY(1:length(find(model_out.trajectoryGlobalFrame.trajectoryPointsX(:,i)))-1,i)).^2).^0.5)];
end
mean_d = mean(distances2(1:end));
max_d = max(distances2(1:end));
min_d = min(distances2(1:end));
std_d = std(distances2(1:end));
disp('-------------------------------------------------')
disp(['Mean of the difference between Sfun and Simulink trajectory ', num2str(mean_d), ' meters']);
disp(['Maximum of the difference between Sfun and Simulink trajectory ', num2str(max_d), ' meters']);
disp(['Minimum of the difference between Sfun and Simulink trajectory ', num2str(min_d), ' meters']);
disp(['Standard deviation of the difference between Sfun and Simulink trajectory ', num2str(std_d), ' meters']);
end

