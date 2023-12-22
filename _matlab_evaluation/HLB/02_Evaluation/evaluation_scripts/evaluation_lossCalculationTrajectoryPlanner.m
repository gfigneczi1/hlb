function [loss] = evaluation_lossCalculationTrajectoryPlanner(traj,reference,GT_U ,GT_Y, ver,P)
if (isempty(traj))
    % there is nothing in the trajectory, we put f = 0
    loss = 0;
else
    ego_path = reference;
    if (ver == 0)
        loss = (abs(ego_path-traj)); 
        loss = sum(sum((loss)));
        loss = loss/(size(ego_path,2));
    elseif (ver == 1)
        loss = 0;
        for i=1:size(traj,1)
            distances = ((ego_path(:,1)-traj(i,1)).^2 + (ego_path(:,2)-traj(i,2)).^2).^0.5;
            idx = find(distances == min(distances),1);
            if (~isempty(idx))
                loss = loss + ((traj(i,1)-ego_path(idx(1,1),1))^2 + (traj(i,2)-ego_path(idx(1,1),2))^2)^0.5;
            end
        end
        loss = loss/size(traj,1);
        %loss = sum(((ego_path(:,1)-traj(:,1)).^2 + (ego_path(:,2)-traj(:,2)).^2).^0.5)/size(traj,1);
    elseif (ver == 2)
        displacement = (diff(traj(:,1)).^2+diff(traj(:,2)).^2).^0.5;
        deviation = ((ego_path(:,1)-traj(:,1)).^2  + (ego_path(:,2)-traj(:,2)).^2).^0.5;
        deviation_corrected = deviation; %.*(exp(U/max(U))-1);
        loss = sum(deviation_corrected)/size(traj,1) + max(deviation)/1.85;
        Ptemp = max(0,sum(P)-1);
        loss = loss + (Ptemp*100);
        Ptemp = min(sum(P)-10/250,0);
        loss = loss + (Ptemp*(-100));
        if ((250-sum(displacement)) > 150)
            loss = loss + 10*(250-sum(displacement))/250;
        end
    elseif (ver==3)
        displacement = (diff(traj(:,1)).^2+diff(traj(:,2)).^2).^0.5;
        deviation = ((ego_path(:,1)-traj(:,1)).^2  + (ego_path(:,2)-traj(:,2)).^2).^0.5;
        deviation_corrected = deviation;
        loss = sum(deviation_corrected)/size(traj,1) + max(deviation)/1.85;
    elseif (ver==4)
        %% Correlation check
        mdl = fitlm(GT_U(1,abs(GT_U(1,:))>0.75),GT_Y(2,abs(GT_U(1,:))>0.75));
        R1 = mdl.Rsquared.Ordinary^0.5;
        mdl = fitlm(GT_U(2,abs(GT_U(2,:))>0.75),GT_Y(3,abs(GT_U(2,:))>0.75));
        R2 = mdl.Rsquared.Ordinary^0.5;
        mdl = fitlm(GT_U(3,abs(GT_U(3,:))>0.75),GT_Y(4,abs(GT_U(3,:))>0.75));
        R3 = mdl.Rsquared.Ordinary^0.5;
        loss = 1/3*(R1+R2+R3);        
    else
        loss = 1000;
    end
end
end


