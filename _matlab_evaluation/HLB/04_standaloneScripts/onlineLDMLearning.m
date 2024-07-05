close all;
clear;

pathData = "C:\git\hlb\_temp\plots\dataOut.mat";
load(pathData);
discPath = "C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\DISC_summary.xlsx";
disc = readtable(discPath);

%% PCA plot
Pcent = [];
for i=1:length(dataOut.Pcentroids)
    Pcent = [Pcent dataOut.Pcentroids{i}];
end
Pdr = [];
for i=1:length(dataOut.data)
    Pdr = [Pdr dataOut.data(i).P_GT];
end
%% DISC plot
% conversion of disc params to +/-1 range X+: B, X-: R, Y+: G, Y-:Y
for i=1:size(disc,1)
    discTransformed(i,:) = [mean([disc.B(i)/25, -disc.R(i)/25]), mean([disc.G(i)/25, -disc.Y(i)/25])];
    plot(discTransformed(i,1), discTransformed(i,2), 'o', 'DisplayName', strcat("Dr", num2str(disc.DrID(i))), 'LineStyle', 'none'); hold on;
end
xlim([-1,1]);
ylim([-1,1]);
xline(0, 'HandleVisibility','off', 'LineWidth',2);
yline(0, 'HandleVisibility','off', 'LineWidth',2);
legend;
grid on;

%% clusters plot
for j=1:length(dataOut.Pcentroids)-1
    for k=j+1:length(dataOut.Pcentroids)
        delta(j,k) = sqrt(sum((dataOut.Pcentroids{j}-dataOut.Pcentroids{k}).^2)/length(dataOut.Pcentroids{k}));
    end
end
delta = delta/max(max(delta));


set(0,'defaulttextInterpreter','latex') ;

config.PARAMETER_CLUSTERING = true;
config.GENERATE_EST_PLOTS = true;
config.USE_SYNTHETIC_DATA = false;
config.DATA_RESTRUCTURE = true;

clc;

for driverId = 1:length(dataOut.data)
    data{driverId} = estimateParameters(driverId, dataOut, config);

    %plotParams(data(driverId).data(1), "left");
    fprintf("Results of right side:\n");
    printResults(data{driverId}.data(1));
    fprintf("\n");
    fprintf("Results of left side:\n");
    printResults(data{driverId}.data(2));
end

if (config.PARAMETER_CLUSTERING)
    driverId = 1;
    for i=1:length(data)
        if (dataOut.data(i).drID ~= 8 && dataOut.data(i).drID ~= 9)
            P(:,driverId) = [reshape(data{driverId}.data(1).GT.Pest,9,1); reshape(data{driverId}.data(2).GT.Pest,9,1); data{driverId}.p0];
            driverId = driverId+1;
        end
    end
    K = 3;
    rng(1); %doc: https://www.mathworks.com/help/matlab/math/controlling-random-number-generation.html
    normalizedParams = P(1:18,:)/max(max(abs(P(1:18,:))));
    normalizedParams = [normalizedParams; P(19,:)/max(abs(P(19,:)))];

    [driverClusters,~,sumD] = kmeans(normalizedParams',K, 'Distance','sqeuclidean', 'MaxIter', 1000, 'Replicates',10); 
    for k=1:K
        Pcentroids.Pcentroids{k} = sum(P(:,driverClusters==k)')'/numel(find(driverClusters==k));
    end
dataOut.Pcentroids = Pcentroids.Pcentroids;
dataOut.driverClusters = driverClusters;
end

classifyDrivers(data, dataOut);

save(fullfile("../../_temp/plots",...
                        'data.mat'), 'data');

%% end of main scripts

%% START of subfunctions

function plotParams(data, token)
    f = figure('Position', [100, 100, 850, 650]);

    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    
    for i=1:size(data.EKF.output,2)
        subplot(3,3,i);
        plot(data.EKF.output(:,i), 'color', 'r', 'DisplayName',strcat("$P_{est,EKF}^{", token, "}$"));
        xlabel('Sample number');
        grid on;
        hold on;
        originalSize = get(gca, 'Position');
        set(gca, 'FontSize', 12);
        if (i==1)
            legend('Location', 'best');            
        end
        set(gca, 'Position', originalSize);
        j = i-floor(i/3);
        k = floor(i/3)+1;

        title(strcat("$P_{", token, "}^{", num2str(k),",", num2str(j),"}$"));

        ylim([max(min(data.EKF.output(:,i))) - 0.5*abs(max(min(data.EKF.output(:,i)))), ...
            min(max(data.EKF.output(:,i))) + 0.5*abs(min(max(data.EKF.output(:,i))))]);
    end


    savefig(f, fullfile("../../_temp/plots",...
                        'Parameter_sweep.fig'));
    saveas(f, fullfile("../../_temp/plots",...
                            'Parameter_sweep.png'));
    close(f);

end

function classifyDrivers(data, dataOut)
% loop through the drivers, do classification and print results
clusters = ["outlier", "dynamic", "balanced"];
driverId = 1;
for i = 1:length(data)
    if (dataOut.data(i).drID ~=8 && dataOut.data(i).drID ~=9)
        P_left = reshape(data{driverId}.data(1).EKF.Pest,3,3);
        P_right = reshape(data{driverId}.data(2).EKF.Pest,3,3);
        P = [reshape(P_left,9,1); reshape(P_right,9,1); data{driverId}.p0];
        %P = [reshape(data{driverId}.data(1).GT.Pest,9,1); reshape(data{driverId}.data(2).GT.Pest,9,1); data{driverId}.p0];
        for clusterId = 1:length(dataOut.Pcentroids)
            Pcent = dataOut.Pcentroids{clusterId};
            deltaP(clusterId) = (sum((Pcent-P).^2)).^0.5;
        end
        [~, c(driverId)] = min(deltaP);
        if (c(driverId)==dataOut.driverClusters(driverId))
            res = "match";
        else
            res = "mismatch";
        end
        fprintf("Cluster for driver %d is: GT: %s, \tEKF: %s, \t res: %s\n", dataOut.data(i).drID, clusters(dataOut.driverClusters(driverId)), clusters(c(driverId)), res);
        driverId = driverId+1;
    end
end
    
end

function printResults(data)
    fprintf("Corresponding parameters are:\n \tR=%f, \n \tQy=%f \n \tQp=%f\n", data.EKF.parameters(1), data.EKF.parameters(2), data.EKF.parameters(3));
    fprintf("Maximum log likelihood: log(L)=%f\n\n", data.EKF.maxLogLikelihood);
    fprintf("Kalman Filter has been initialized with the following values:\n");
    fprintf("\t x0 = [%f %f %f %f]\n", data.EKF.initValues.x0);
    fprintf("\t P0 = diag([%f %f %f %f]\n", diag(data.EKF.initValues.P0));
    fprintf("Final NRMS values:\n");
    fprintf("\tEKF:\t%f\n", data.EKF.terminalError);
    MA_conv = convergenceCheck(data.MA.eL2, 0.002, 50);
    EKF_conv = convergenceCheck(data.EKF.eL2, 0.002, 50);
    fprintf("Convergence is reached:\n");
    fprintf("\tEKF: %d, relative to MA: %f\n", EKF_conv, EKF_conv/MA_conv);
    fprintf("\n");
    fprintf("Run time:\n");
    fprintf("\tEKF: %f secs, relative to MA: %f\n", mean(data.EKF.dt), mean(data.EKF.dt)/mean(data.MA.dt));
end

function [convIndex] = convergenceCheck(data, P_conv, P_minConvLength)
% convergence is calculated based on the gradient of the data
% input parameter: absolute change in value
% if convergence is reached and keeps up till the end, the convergence is
% true. If convergence is reached for only a small amount of time, it is
% considered to be non-converged.
gradients = movmean(diff(movmean(data,100)),100);
if (isempty(find(abs(gradients(end:-1:1)) >= P_conv,1)))
    convIndex = inf;
else
    convIndex = length(data) - find(abs(gradients(end:-1:1)) >= P_conv,1);
    if ((length(data)-convIndex) <= P_minConvLength)
        convIndex = inf;
    end
end
end

function driverData = estimateParameters(driverID, dataOut, config)
    %% In the followings
    % Several methods are compared to each other to estimate the parameters of
    % LDM. Estimation mainly refers to the P matrices (curve matrices)
    % estimation. First, the straight line offset parameter is calculated, then
    % the straight line offset is subtracted from the overall data to provide
    % proper inputs for curve parameter estimation.

    driverData.driver = dataOut.data(driverID).drID;
    
    % regression happens for first raw of parameters
    U = dataOut.data(driverID).U;
    dY = dataOut.data(driverID).dY;

    [P_GT_, U_left, dY_left, U_right, dY_right] = functional_driverModelLearning(dataOut.data(driverID).U', dataOut.data(driverID).dY',8);

    P_GT = dataOut.data(driverID).P_GT;
    driverData.p0 = P_GT(19);
    
    %% solution 1: linear regression
    
    % restructure data
    if (config.DATA_RESTRUCTURE)
        [U_left, dY_left] = dataRestructure (U_left, dY_left, "random");
        [U_right, dY_right] = dataRestructure (U_right, dY_right, "random");
    end
    
    tic;
    P_LDM_left = reshape(P_GT(1:9), 3,3);
    
    driverData.data(1).MA.dt = toc;
    driverData.data(2).MA.dt = driverData.data(1).MA.dt;

    P_LDM_right = reshape(P_GT(10:18), 3,3);
    clear P_GT;

    for side=1:2
        global U dY P_posterior0 P_LDM
        if (side==1)
            P_LDM = P_LDM_left;
            U = U_left;
            dY = dY_left;
            token = 'left';            
        else
            P_LDM = P_LDM_right;
            U = U_right;
            dY = dY_right;
            token = 'right';
        end
        
        % EKF estimator
        % optimization

        P_posterior0 = eye(12)*10;
        [xCons,fCons] = fmincon(@maxLogLikelihood, [0, 0, 0], [],[],[],[],zeros(1,3),[inf inf inf], []);
                
        sy = xCons(1);
        sp = xCons(2);
        ry = xCons(3);
        Q = [eye(3)*sy^2 zeros(3,9); zeros(9,3) eye(9)*sp^2]; 
    
        % Kalman-filter initialization
        Pinit = zeros(9,1); %estimateInitialParameters(U, dY);
        R = eye(3)*ry^2; 
    
        H = zeros(3,12);
        H(:,1:3) = eye(3);
    
        x_posterior = [dY(1,:)'; Pinit];
        P_posterior = P_posterior0;
    
        x0 = x_posterior;
        
        sigma_xx = P_posterior0;
        sigma_yy = H*sigma_xx*H'+R;
    
        Yerror = dY(1,:)'-H*x_posterior;
        L = log(det(sigma_yy))+Yerror'*inv(sigma_yy)*Yerror;
    
        for i=1:size(U,1) 

            tic;
            G = [x_posterior(4:6)'; x_posterior(7:9)'; x_posterior(10:12)'; zeros(9,3)];
            F = [zeros(3,3) [U(i,:) zeros(1,6); zeros(1,3) U(i,:) zeros(1,3); zeros(1,6) U(i,:)]; [zeros(9,3) eye(9)]];
    
            [x_posterior, P_posterior] = ekf_simple(x_posterior, P_posterior, U(i,:)', F,H,dY(i,:)',R,Q);

            parameters(i,1:9) = x_posterior(4:12);
    
            sigma_xx = F*sigma_xx*F'+Q;
            sigma_yy = H*sigma_xx*H'+R;

            dtEkf(i) = toc;
    
            Yerror = dY(i,:)'-H*x_posterior;
            L = L+log(det(sigma_yy))+Yerror'*inv(sigma_yy)*Yerror;
                
            P_GT = functional_driverModelLearning(U(1:i,:), dY(1:i,:), 8); %inv(U(1:i,:)'*U(1:i,:))*U(1:i,:)'*dY(1:i,:); 
            if (side==1)
                P_GT = reshape(P_GT(1:9),3,3);
            else
                P_GT = reshape(P_GT(10:18),3,3);
            end
    
            P_est = [x_posterior(4:6) x_posterior(7:9) x_posterior(10:12)];
    
            normP = (max(max(P_LDM))-min(min(P_LDM)));
            eL2_post_LDM(i) = (sum(sum((P_est-P_LDM).^2))/9)^0.5 / normP;
            eL2_moving_LDM(i) = (sum(sum(((P_GT-P_LDM)).^2))/9)^0.5 / normP;
        end %end of KF estimation via data

        L = 1/2*L;
    
        if (config.GENERATE_EST_PLOTS)
            f = figure('Position', [100, 100, 750, 420]);
            set(f,'defaulttextInterpreter','latex') ;
            set(f, 'defaultAxesTickLabelInterpreter','latex');  
            set(f, 'defaultLegendInterpreter','latex');
            subplot(2,1,1);                
            plot(reshape(P_LDM,9,1), 'color', 'k', 'Marker', 'x', 'DisplayName','$P_{GT}^{x,1}$', 'MarkerSize', 10);
            hold on;
            plot(reshape(P_GT,9,1), 'color', 'b', 'Marker', 'o', 'DisplayName','$P_{est,MA}^{x,1}$', 'MarkerSize', 10);
            plot(x_posterior(4:12), 'color', 'r', 'Marker', 'x', 'DisplayName','$P_{est,EKF}^{x,1}$', 'MarkerSize', 10);
            grid on;
            xticks([0,1,2,3,4,5,6,7,8,9,10]);
            xticklabels(["", strcat("$P_{",token,"}^{1,1}$"), strcat("$P_{",token,"}^{2,1}$"), strcat("$P_{",token,"}^{3,1}$"), ...
                strcat("$P_{",token,"}^{1,2}$"), strcat("$P_{",token,"}^{2,2}$"), strcat("$P_{",token,"}^{2,3}$"), ...
                strcat("$P_{",token,"}^{1,3}$"), strcat("$P_{",token,"}^{2,3}$"), strcat("$P_{",token,"}^{3,3}$"), ""]);
            xlim([0,10]);
            set(gca,'FontSize',12);
            title(strcat('\textbf{Final Parameter distribution - $P_{', token, '}$}'), 'FontSize', 14);
    
            legend('Location', 'best', 'FontSize', 12);
    
            subplot(2,1,2);
            plot(movmean(eL2_moving_LDM,20), "DisplayName", 'NRMS($\Delta P_{MA}^{GT}$)', 'color', 'b', 'LineWidth',2); hold on;
            plot(movmean(eL2_post_LDM,20), "DisplayName", 'NRMS($\Delta P_{est,EKF}^{GT}$)', 'color', 'r', 'LineWidth',2);
    
            grid on;
            title(strcat('\textbf{NRMS value of error between parameters - $P_{', token, '}$}'), 'FontSize', 14);
            xlabel('Sample number', 'FontSize', 12);
            ylabel('NRMS($\Delta P$)', 'FontSize', 12);
            legend('Location', 'best', 'FontSize', 12);
            ylim([0,2]);
            yline(1, 'LineWidth',2, 'color', [150 150 150]/255, 'HandleVisibility','off');
            yline(0.15, 'LineWidth',1, 'LineStyle', '--', 'color', 'k', 'label', "$conv_{thd}=15\%$", 'Interpreter', 'latex', 'HandleVisibility','off');
    
            set(gca,'FontSize',12);
            xlim([0, length(eL2_post_LDM)]);
    
            savefig(f, fullfile("../../_temp/plots",...
                        strcat('Estimation_', 'Driver_', num2str(driverID),'_side_', token, '.fig')));
            saveas(f, fullfile("../../_temp/plots",...
                            strcat('Estimation_', 'Driver_', num2str(driverID),'_side_', token, '.png')));
            close(f);
        end  
        driverData.data(side).EKF.maxLogLikelihood = L;
        driverData.data(side).EKF.terminalError = eL2_post_LDM(end);
        driverData.data(side).EKF.parameters = [ry, sy, sp];
        driverData.data(side).EKF.initValues.x0 = x0;
        driverData.data(side).EKF.initValues.P0 = P_posterior0;
        driverData.data(side).EKF.eL2 = eL2_post_LDM;
        driverData.data(side).EKF.Pest = x_posterior(4:12);
        driverData.data(side).MA.eL2 = eL2_moving_LDM;
        driverData.data(side).MA.Pest = P_GT;
        driverData.data(side).GT.Pest = P_LDM;
        driverData.data(side).EKF.output = parameters;
        driverData.data(side).EKF.dt = dtEkf;

        clear eL2_moving_LDM eL2_post_LDM NLMS LMS nrms nrmsLms parameters parameters_nlms dtEkf dtLms dtNlms
    end % end of side separation
end

function f = maxLogLikelihood(x0)
global U dY P_posterior0 P_LDM
    % initialize state
    sy = x0(1);
    sp = x0(2);
    ry = x0(3);
    Q = [eye(3)*sy^2 zeros(3,9); zeros(9,3) eye(9)*sp^2];

    % Kalman-filter initialization
    Pinit = zeros(9,1); %estimateInitialParameters(U, dY);
    R = eye(3)*ry^2;

    H = zeros(3,12);
    H(:,1:3) = eye(3);

    x_posterior = [dY(1,:)'; Pinit];
    sigma_xx = P_posterior0;
    sigma_yy = H*sigma_xx*H'+R;
    P_posterior = P_posterior0;

    Yerror = dY(1,:)'-H*x_posterior;
    L = log(det(sigma_yy))+Yerror'*inv(sigma_yy)*Yerror;
    eL2_post_LDM = 0;
   
    for i=1:size(U,1)    
        G = [x_posterior(4:6)'; x_posterior(7:9)'; x_posterior(10:12)'; zeros(9,3)];
        F = [zeros(3,3) [U(i,:) zeros(1,6); zeros(1,3) U(i,:) zeros(1,3); zeros(1,6) U(i,:)]; [zeros(9,3) eye(9)]];

        [x_posterior, P_posterior] = ekf_simple(x_posterior, P_posterior, U(i,:)', F,H,dY(i,:)',R,Q);

        sigma_xx = F*sigma_xx*F'+Q;
        sigma_yy = H*sigma_xx*H'+R;

        Yerror = dY(i,:)'-H*x_posterior;
        L = L+log(det(sigma_yy))+Yerror'*inv(sigma_yy)*Yerror;

        P_est = [x_posterior(4:6) x_posterior(7:9) x_posterior(10:12)];
    
        normP = (max(max(P_LDM))-min(min(P_LDM)));
        eL2_post_LDM = eL2_post_LDM+(sum(sum((P_est-P_LDM).^2))/9)^0.5 / normP;
    end
    L = -1/2*L;
    f = -L; %eL2_post_LDM/size(U,1);
end

function P0 = estimateInitialParameters(U, dY)
U00 = U(1:40,:);
dY00 = dY(1:40,:);

% P0 is calculated from a subset of measured data
P0 = inv(U00'*U00)*U00'*dY00;
P0 = reshape(P0,9,1);

end

function [y, Cyy, Cxy]= UT(H, x, P, R)

N = size(x,1);
for i=1:2*N+1
    if (i==1)
        X(1:N,i) = mean(x);
    elseif (i<=N+1)
        Cxx = (N*P).^0.5;
        X(:,i) = mean(x)+Cxx(:,i-1);
    else
        Cxx = (N*P).^0.5;
        X(:,i) = mean(x)-Cxx(:,i-N-1);
    end
end

Z = H*X;
W = [zeros(size(H,1),1) ones(size(H,1),2*N)*1/(2*N)];

y = 0; Cyy = R; Cxy = 0;

for i=1:2*N+1
    y = y+W(i)*Z(:,i);
end

for i=1:2*N+1
    Cyy = Cyy + W(i)*(Z(:,i)-y)*(Z(:,i)-y)';
    Cxy = Cxy +W(i)*(X(:,i)-mean(x))*(Z(:,i)-y)';
end

end

function [U, dY] = dataRestructure (U, dY, mode)

if (mode=="ascend" || mode=="descend")
     % adding extra as the mean
    U = [U mean(U,2)];
    dataSorted = sortrows([U dY],4, mode);
    dataShuffled = zeros(size(dataSorted));
    dataShuffled(1:2:2*round(size(dataSorted,1)/2),:) = dataSorted(1:round(size(dataSorted,1)/2),:);
    dataShuffled(2:2:round(size(dataSorted,1)),:) = dataSorted(round(size(dataSorted,1)/2)+1:end,:);
elseif (mode=="random")
    shuffledIndeces = randperm(size(U,1));
    dataShuffled(:,1:3) = U(shuffledIndeces,:);
    dataShuffled(:,5:7) = dY(shuffledIndeces,:);
end

U = dataShuffled(:,1:3);
dY = dataShuffled(:,5:7);
end

function [U, dY] = createSyntheticData (U, dY, noiseMode, noiseFactor)
% this function calculates some P matrices, then produces ground truth
% data, and also adds noise per the inputs 'noiseMode' and 'noiseFactor'
P = inv(U'*U)*U'*dY(:,:);
dY = U*P;
if (noiseMode == "randomProportional")
    noise = rand(size(dY))*2-1;
    normalizedNoiseFactor = noiseFactor*abs(dY);
    dY = dY + noise.*normalizedNoiseFactor;
end

end

function [x_posterior, P_posterior] = ekf_simple(x_posterior, P_posterior, u, F,H,z,R,Q)
    % prediction
    x_prior = stateUpdate(x_posterior, u');
    P_prior = F*P_posterior*F'+Q;

    % correction
    y = outputCalculation(x_prior);
    y_error = z-y;
    K = P_prior*H'*inv(H*P_prior*H'+R);
    x_posterior = x_prior+K*y_error;
    P_posterior = (eye(12)-K*H)*P_prior;

end

function x = stateUpdate(x,u)
    P = reshape(x(4:12), 3,3);
    x(1:3) = u*P;
end

function y = outputCalculation(x)
    y = x(1:3);
end