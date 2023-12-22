clear; close all;

input_folder = "../../_temp";

% load mat files
matFiles = dir(fullfile(input_folder,'/*.mat'));

for fileID = 1:length(matFiles)
    rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
    drID = matFiles(fileID).name(1:5);
    
    [x,y,utmzone] = deg2utm(rawData.LatPos_abs,rawData.LongPos_abs);
    theta_calc = zeros(length(x),1);
    theta_calc(1) = 0;
    for i=2:length(x)
        if (x(i) == x(i-1))
            if (y(i) == y(i-1))
                theta_calc(i) = theta_calc(i-1);
            elseif (y(i) > y(i-1))
                theta_calc(i) = pi()/2;
            else
                theta_calc(i) = -pi()/2;
            end
        elseif (y(i) == y(i-1))
            if (x(i) > x(i-1))
                theta_calc(i) = 0;
            else
                theta_calc(i) = -pi();
            end
        else
            theta_calc(i) = atan((y(i)-y(i-1))/(x(i)-x(i-1)));
            if (y(i) < y(i-1))
                if (x(i) < x(i-1))
                    theta_calc(i) = theta_calc(i) + pi();
                end 
            else
                if (x(i) < x(i-1))
                    theta_calc(i) = theta_calc(i) + pi();
                end
            end                       
        end
    end

    % creating the corridor information
    for j=1:size(x,1)
        map(j,1:2) = pos_tf2GPS(x(j),y(j),theta_calc(j),rawData.c01_left(j));
        map(j,3:4) = pos_tf2GPS(x(j),y(j),theta_calc(j),rawData.c01_right(j));
    end
    
    mapFixed = map((rawData.c01_left~=0 .* rawData.c01_right~=0), :);
    csvwrite(fullfile(matFiles(fileID).folder,strcat(matFiles(fileID).name(1:end-4), 'map.csv')),mapFixed);
end