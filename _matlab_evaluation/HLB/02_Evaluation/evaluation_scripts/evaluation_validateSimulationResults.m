function evaluation_validateSimulationResults(traj, output)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    hold on
    plot(traj(:,1),traj(:,2));
    
    modelTrajX = squeeze(output.trajectoryGlobalFrame(1,1,:));
    modelTrajY = squeeze(output.trajectoryGlobalFrame(1,2,:));
    plot(modelTrajX(:,1),modelTrajY(:,1))
    
    k = 0;
    for i=1:length(modelTrajX)-1
        if (modelTrajX(i-k) == modelTrajX(i-k+1)) && (modelTrajY(i-k) == modelTrajY(i-k+1))
            modelTrajX(i-k) = [];
            modelTrajY(i-k) = [];
            k = k + 1;
        end
    end
    
    modelTrajY_alligned = interp1(modelTrajX, modelTrajY, traj(:,1), 'spline', 'extrap');
   
    error = 0;
    for i=1:length(traj(:,1))
        error = error + abs(modelTrajY_alligned(i) - traj(i,2));
    end
    errorPerStep = error/length(traj(:,2));
   
end

