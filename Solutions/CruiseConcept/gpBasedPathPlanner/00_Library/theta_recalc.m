function [theta_calc] = theta_recalc(segment, correction)
%THETA Summary of this function goes here
%   Detailed explanation goes here
    X_abs = movmean(segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360, 1);
	Y_abs = movmean(segment.LatPos_abs * 111.32*1000,1);

    % recalculating theta base
%     theta_calc = zeros(length(X_abs),1);
%     theta_calc(1) = 0;
%     for i=2:length(X_abs)
%         if (X_abs(i) == X_abs(i-1))
%             if (Y_abs(i) == Y_abs(i-1))
%                 theta_calc(i) = theta_calc(i-1);
%             elseif (Y_abs(i) > Y_abs(i-1))
%                 theta_calc(i) = pi()/2;
%             else
%                 theta_calc(i) = -pi()/2;
%             end
%         elseif (Y_abs(i) == Y_abs(i-1))
%             if (X_abs(i) > X_abs(i-1))
%                 theta_calc(i) = 0;
%             else
%                 theta_calc(i) = -pi();
%             end
%         else
%             theta_calc(i) = atan((Y_abs(i)-Y_abs(i-1))/(X_abs(i)-X_abs(i-1)));
%             if (Y_abs(i) < Y_abs(i-1))
%                 if (X_abs(i) < X_abs(i-1))
%                     theta_calc(i) = theta_calc(i) + pi();
%                 end 
%             else
%                 if (X_abs(i) < X_abs(i-1))
%                     theta_calc(i) = theta_calc(i) + pi();
%                 end
%             end                       
%         end
%     end
    
if (correction == 1)
    theta_calc = diff(segment.Y_abs)./diff(segment.X_abs);
    if (size(theta_calc,1) == 1)
        theta_calc = [theta_calc(1) theta_calc];
    else
        theta_calc = [theta_calc(1); theta_calc];
    end
    theta_calc = atan(theta_calc); %+pi();
else
    theta_calc = diff(segment.Y_abs)./diff(segment.X_abs);
    if (size(theta_calc,1) == 1)
        theta_calc = [theta_calc(1) theta_calc];
    else
        theta_calc = [theta_calc(1); theta_calc];
    end
    theta_calc = atan(theta_calc)+pi();
end

end

