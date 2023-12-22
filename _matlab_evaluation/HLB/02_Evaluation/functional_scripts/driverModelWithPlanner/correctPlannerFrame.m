function plannerFrameCorrected = correctPlannerFrame(X_abs, Y_abs, plannerFrame)
plannerFrameCorrected = plannerFrame;
% searching for closest point for X-Y
distances = ((X_abs-plannerFrame(1)).^2+(Y_abs-plannerFrame(2)).^2).^0.5;
idx = find(distances == min(distances),1);

if(~isempty(idx))
    if (distances(idx(1,1)) >=1.25)
        plannerFrame(1) = X_abs(idx(1,1))-sin(plannerFrame(3))*1.25;
        plannerFrame(2) = Y_abs(idx(1,1))+cos(plannerFrame(3))*1.25;
    end
end        

end

