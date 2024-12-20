function [estimation, npDistances, input] = generateEstimate(modelID,inputRaw)

global parameters params

%% step 1: norm and central
for i=1:size(inputRaw,2)
    input(:,i) = (inputRaw(:,i)-params.c_in(i))/params.s_in(i);
end

%% step 2: model based estimate generation
switch modelID
    case "sgp"
        meanfunc = [];       % Start with a zero mean prior
        eval(strcat('covfunc = ',parameters.PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
        % ID problem            
        likfunc = @likGauss;    % Gaussian likelihood
        for i = 1:length(params.hyp_opt_array)
            hyp = struct('mean', [], 'cov', 0, 'lik', -1);
            hyp.cov = params.hyp_opt_array{i}.cov;
            hyp.lik = params.hyp_opt_array{i}.lik;
            [estimationNormed(:,i), deviationNormed(:,i)] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, params.input_estimation{i}, params.output_estimation{i}, input); % extract the mean and covarance functions
        end

        npDistances = parameters.PARAMS.OUTPUT_SHIFT;
    case "ldm"
        for n=1:size(params.P,1)
            estimationNormed(:,n) = input*params.P(n,:)';
        end
        npDistances = parameters.PARAMS.npDistances;
    case "eldm"
        for n=1:size(params.P_left,1)
            estimationNormed(:,n) = input*[params.P_left(n,:)'; params.P_right(n,:)'];
        end
        npDistances = parameters.PARAMS.npDistances;
    case "lrm"
        for n=1:size(params.P,1)
            estimationNormed(:,n) = input*params.P(n,:)';
        end
        npDistances = parameters.PARAMS.OUTPUT_SHIFT;
    case "phtpm"
        inputCell{1,1} = input';
        estimationNormed = predict(params.net, inputCell);
        estimationNormed = double(estimationNormed{1,1});
        dx = inputRaw(end,3)*0.1;
        npDistances = dx:dx:size(estimationNormed,2)*dx;
        estimationNormed = estimationNormed(end:-1:1);
        
end

%% step 3: estimation denormalization and decentralization
for i=1:size(estimationNormed,2)
    estimation(:,i) = estimationNormed(:,i)*params.s_out(i);
    if (~parameters.eliminateOffset)
        estimation(:,i) = estimation(:,i)+params.c_out(i);
    end
end

end

