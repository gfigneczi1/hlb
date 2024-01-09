function [transformed_coefficients] = lane_border_interpolation_2ndorder(base,client,coefficients, dtheta_in, variation)
%LANE_BORDER_INTERPOLATION_2NDORDER Summary of this function goes here
%   Base: [X Y theta] in absolute frame of base coordinate system origo
%   (same for client) - NOTE: client is ahead of base in time
%   Coefficients: [c0 c1 c2] in base coordinate system
    dtheta = dtheta_in;
    % Due to max 2pi angle limitation on GPS side, a slight correction of
    % dtheta is necessary, whenever base turns exactly 2pi radian compared
    % to client. If that is the case, dtheta must be equal to 0 radian, as
    % the angle was continuous and the 2pi radian turn within 80ms (video
    % cycle time) is absolutely not feasible.
    eps = 1 * pi()/180; % tolerance is 1 degree
    if ((abs(dtheta) <= eps/2) && (abs(dtheta) >= -eps/2))
        v2 = inf;
    elseif dtheta == pi()
        v2 = -inf;
    else
        v2 = tan(dtheta - pi()/2);
        v2 = v2/(v2^2+1^2)^0.5;
        v1 = (1-v2^2)^0.5;
        v = v2/v1;
    end
    dx = client(1) - base(1);
    dy = client(2) - base(2);
    
    x0 = cos(base(3))*dx + sin(base(3))*dy;
    y0 = -sin(base(3))*dx + cos(base(3))*dy;
    
    if (v2 == inf)
        % coordinate systems are parallel
        x1 = x0;
        y1 = coefficients(1)+coefficients(2)*x1+coefficients(3)*x1^2 + y0;
        c0v = y1;
        c1v = coefficients(2)+2*coefficients(3)*x1; % orientation does not change only due to translation
        c2v = coefficients(3); % curvature does not change
    elseif (v2 == -inf)
        % coordinate systems are parallel, but pointing against each other
        x1 = x0;
        y1 = coefficients(1)+coefficients(2)*x1+coefficients(3)*x1^2 + y0;
        c0v = y1;
        c1v = coefficients(2)+2*coefficients(3)*x1; % orientation does not change only due to translation
        c2v = coefficients(3); % curvature does not change
        c0v = -c0v;
        c2v = -c2v;
    else
        % coordinate systems are not parallel
        % calculating POI as the crossing point of the line and the client
        % Y axis
        D = (coefficients(2)-v)^2-4*coefficients(3)*(coefficients(1)+v*x0-y0);
        if (D<0)
            c0v = 0; c1v = 0; c2v=0;
        elseif (D==0)
            x1 = -(coefficients(2)-v)/(2*coefficients(3));
            y1 = coefficients(1)+coefficients(2)*x1+coefficients(3)*x1^2;
            % now x1 and y1 are in the base frame, now transform to client
            % frame
            x1base = cos(dtheta)*x1 + sin(dtheta)*y1;
            y1base = -sin(dtheta)*x1 + cos(dtheta)*y1 - y0;
            c0v = y1base;
            c1v = coefficients(2)+2*coefficients(3)*x1 - dtheta;
            c2v = coefficients(3); % curvature does not change
        else
            x11 = (-(coefficients(2)-v)+((coefficients(2)-v)^2-4*coefficients(3)*(coefficients(1)+v*x0-y0))^0.5)/(2*coefficients(3));
            x12 = (-(coefficients(2)-v)-((coefficients(2)-v)^2-4*coefficients(3)*(coefficients(1)+v*x0-y0))^0.5)/(2*coefficients(3));
            if (abs(x11) <= abs(x12))
                x1 = x11;
            else
                x1 = x12;
            end
            y1 = coefficients(1)+coefficients(2)*x1+coefficients(3)*x1^2;
            % now x1 and y1 are in the base frame, now transform to client
            % frame
            x1base = cos(dtheta)*x1 + sin(dtheta)*y1;
            y1base = -sin(dtheta)*x1 + cos(dtheta)*y1 - y0;
            c0v = y1base;
            c1v = coefficients(2)+2*coefficients(3)*x1 - dtheta;
            c2v = coefficients(3); % curvature does not change
        end
    end
    transformed_coefficients = [c0v c1v c2v];
end

