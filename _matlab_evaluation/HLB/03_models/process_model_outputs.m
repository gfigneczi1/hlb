input = 0;
switch input
    case 0
        clear tGX tGY cGX cGY
        tGX(:) = out.model_out.trajectoryGlobalFrame.Data(1,1,:);
        tGY(:) = out.model_out.trajectoryGlobalFrame.Data(1,2,:);
        cGX(:) = out.model_out.corridorGlobalFrame.corridorGlobal_X.Data;
        cGY(:) = out.model_out.corridorGlobalFrame.corridorGlobal_Y.Data;
        % local plots - planner frame
        for i=40:length(out.model_out.corridorPlannerFrame.corridor_X.Time)
            n  = out.model_out.corridorPlannerFrame.corridor_X.Data(i,301);
            cLX = out.model_out.corridorPlannerFrame.corridor_X.Data(i,1:n);
            cLY = out.model_out.corridorPlannerFrame.corridor_Y.Data(i,1:n);
            n  = out.model_out.trajectoryPlannerFrame.Data(3001,1,i);
            tLX = out.model_out.trajectoryPlannerFrame.Data(1:n,1,i);
            tLY = out.model_out.trajectoryPlannerFrame.Data(1:n,2,i);
            X = out.model_out.segments.Data(1:3,1,i);
            Y = out.model_out.segments.Data(1:3,2,i);
            subplot (3,1,1);
            Xego = out.model_out.egoPosePlannerFrame.Data(i,1);
            Yego = out.model_out.egoPosePlannerFrame.Data(i,2);
            thetaego = out.model_out.egoPosePlannerFrame.Data(i,3);
            Xinit = out.model_out.egoPosePlannerFrameInit.Data(i,1);
            Yinit = out.model_out.egoPosePlannerFrameInit.Data(i,2);
            thetainit = out.model_out.egoPosePlannerFrameInit.Data(i,3);
            % ego to planner frame
            Xplanner = (cos(thetainit)*Xego - sin(thetainit)*Yego)+Xinit;
            Yplanner = (sin(thetainit)*Xego + cos(thetainit)*Yego)+Yinit;
            thetaplanner = thetainit + thetaego;
            x = linspace(0,4,10)*cos(thetaplanner);
            y = linspace(0,4,10)*sin(thetaplanner);
            
            plot(cLX,cLY,'bo',tLX,tLY,'kx',Xplanner,Yplanner,'ko',Xplanner+x,Yplanner+y,X,Y,'rx');
            grid on;
            title(num2str(out.model_out.egoPosePlannerFrameInit.Data(i,:)));
            
            
            n  = out.model_out.corridorEgoFrame.corridorEgo_X.Data(i,301);
            cLX = out.model_out.corridorEgoFrame.corridorEgo_X.Data(i,1:n);
            cLY = out.model_out.corridorEgoFrame.corridorEgo_Y.Data(i,1:n);
            n  = out.model_out.trajectoryEgoFrame.Data(3001,1,i);
            tLX = out.model_out.trajectoryEgoFrame.Data(1:n,1,i);
            tLY = out.model_out.trajectoryEgoFrame.Data(1:n,2,i);
            
            subplot(3,1,2);
            plot(cLX,cLY,'bo',tLX,tLY,'kx');
            grid on;
            title(num2str(out.model_out.trajectoryEgoFrame.Time(i)));
            
            subplot(3,1,3);
            
            pause(0.1);
        end
end