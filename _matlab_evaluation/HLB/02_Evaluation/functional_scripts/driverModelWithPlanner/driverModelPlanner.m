function [X, Y, theta, indeces, nominal, U] = driverModelPlanner(segment, P, corridor, P3in, ver)
% This function is the Linear Driver Model which calculates the node point
% offsets based on the input corridor information. The model is
% parametrized through the P matrix. 
% Inputs:
% - segment: subsegment which is cut from the global segment based on the
% ego pose global, containing data velocityX, yawRate. Currently this is
% only a placeholder and not used within the code. In the future this might
% by used, therefore we keep it as an input.
% - corridor: n*4 array, columns are: corridor_X, corridor_Y (midlane
% points 2D coordinates in the planner frame, unit: m), corridor orientation (in the
% planner frame, unit: radian) and corridor curvature (unit: 1/m^2).
% - P: 21*1 parameter vector, which is arranged to 3/7 matrix in this
% function. These are the linear coefficients of the LDM.
% - P3in: nominal road point distances.
% Outputs:
% - X,Y: node point coordinates in the 2D planner frame (4 node points, 1st one
% is always 0,0. (unit: m)
% - theta: node point orientation in the planner frame, (unit: radian)
% - indeces: indeces of the nodepoints based on the corridor points
% (between 1 and the number of the corridor points)
% - nominal: nominal road points extracted from the corridor, in the
% planner frame.
% - U: input vector, which is 7*1 vector, 1-3: average curvature, 4-7:
% curvature gradients in the node point distances.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-07

global GT_U_GP alpha Kxx

%% Introduction
% reshape the parameter matrix

Np = length(P3in); % number of node points is calculated based on the distances given in this vector

    %% Nominal road point calculation
    % first determine the sections within the subsegment. This will yield the
    % nominal (corridor) coordinates in the frame which it receives the
    % corridor info.
    [indeces, ~] = nodePointModel(corridor, P3in);
    X_nominal = corridor(indeces,1);
    Y_nominal = corridor(indeces,2);
    theta_nominal = corridor(indeces,3);

    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(corridor(indeces(j):indeces(j+1),4));
    end

    dkappa_vector = diff(corridor(:,4))./diff(corridor(:,1));
    dkappa_vector = [dkappa_vector(1); dkappa_vector];
    nominal = [X_nominal; Y_nominal; theta_nominal];
    dkappa_nominal = dkappa_vector(indeces);
    
        %% Model calculation
        if (ver==0)
            % LDM
            % initialize variables
            X = zeros(Np+1,1);
            Y = zeros(Np+1,1);
            theta = zeros(Np+1,1);
            % parameters
            P1 = reshape(P(1:Np*Np,1),Np, Np);
            
            U = [kappa_nominal];

            % scaling of the predictor variables
            kappa_lim = 0.001;
            dkappa_lim = 1e-5;
            U_lim = [ones(Np,1)*kappa_lim];
            U = U./U_lim;

            % The ouput of the model is the linear combination of predictor variables
            x = P1*U;
            x(1:Np) = min(max(-1.25,x(1:Np)), 1.25);
        elseif (ver == 1)
            % E-LDM
            kappa_lim = 0.001;
            U = [max(kappa_nominal,0); min(kappa_nominal,0); 1];
            U_lim = [ones(3,1)*kappa_lim; ones(3,1)*kappa_lim; 1];
            U = U./U_lim;
            x = P1*U;
            x(1:3) = min(max(-1.25,x(1:3)), 1.25);
        elseif (ver == 2)
            global GP_mode
            % 0: mean, 1: upper variance, 2: lower variance
            % Gaussian process
            kappa_lim = 0.001;
            U = [kappa_nominal; dkappa_nominal];
            U_lim = [ones(3,1)*kappa_lim; ones(3,1)*kappa_lim; 1];
            U = U./U_lim;
            
            
            sigma_ker=0.2;
            sigma_es = 0.15;
            kernel = @(x1,x2)  exp(-(x1-x2)'*(x1-x2)/(2*sigma_ker));
            x_test=kappa_nominal/kappa_lim;
            N_test=size(x_test,2);

            % Standard implementation
            Ktest=[];
            for i=1:N_test
                for j=1:size(GT_U_GP,2)
                    Ktest(i,j)=kernel(x_test(:,i),GT_U_GP(1:3,j));
                end
            end
            f_est=Ktest*alpha;
            
            
            if (GP_mode == 0)
                x(1:3) = min(max(-1.25,f_est'), 1.25);
                x = x';
            elseif (GP_mode == 1)
                for i=1:N_test
                    m=(Kxx+sigma_es*eye(size(GT_U_GP,2)))\Ktest(i,:)';
                    m=Ktest(i,:)*m;
                    m=m'*m;

                    f_var(i)=kernel(x_test(i),x_test(i))- m; 
                end
                x(1:3) = min(max(-1.25, 2*sqrt(f_var)+f_est),1.25);
                x = x';
            elseif (GP_mode == -1)
                for i=1:N_test
                    m=(Kxx+sigma_es*eye(size(GT_U_GP,2)))\Ktest(i,:)';
                    m=Ktest(i,:)*m;
                    m=m'*m;

                    f_var(i)=kernel(x_test(i),x_test(i))- m; 
                end
                x(1:3) = min(max(-1.25, -2*sqrt(f_var)+f_est),1.25);
                x = x';
            end
        elseif (ver==3)
            global hyp_opt1 hyp_opt2 hyp_opt3 GT_Y_GP
            kappa_lim = 0.001;
            x_test=kappa_nominal/kappa_lim;
            meanfunc = [];          % Start with a zero mean prior
            covfunc = @covSEiso;    % Squared Exponental covariance function
            %% ID problem
            likfunc = @likGauss;    % Gaussian likelihood
            [m1,~] = gp(hyp_opt1, @infGaussLik, meanfunc, covfunc, likfunc, GT_U_GP(1:3,:)', GT_Y_GP(2,:)', x_test'); % extract the mean and covarance functions
            [m2,~] = gp(hyp_opt2, @infGaussLik, meanfunc, covfunc, likfunc, GT_U_GP(1:3,:)', GT_Y_GP(3,:)', x_test'); % extract the mean and covarance functions
            [m3,~] = gp(hyp_opt3, @infGaussLik, meanfunc, covfunc, likfunc, GT_U_GP(1:3,:)', GT_Y_GP(4,:)', x_test'); % extract the mean and covarance functions
            x(1:3) = min(max(-1.25, [m1 m2 m3]),1.25);
            x = x';
        end
        Y(2:Np+1) = Y_nominal(2:Np+1) + x(1:Np); %cos(theta_nominal(2:4)).*x(1:3);

% The outputs in the local coordinate frame are initialized in [0 0 0],
% therefore their first index is not assigned to anything but kept at zero.
X(2:Np+1) = X_nominal(2:Np+1); % - sin(theta_nominal(2:4)).*x(1:3);
theta(2:Np+1) = atan(tan(theta_nominal(2:Np+1))); % limiting the orientation
end

