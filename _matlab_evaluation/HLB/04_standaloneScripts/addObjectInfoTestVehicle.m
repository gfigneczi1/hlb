clear;
input_folder =  'C:\git\kdp_hlb_evalframework\_temp';
output_folder = "C:\git\kdp_hlb_evalframework\_temp";

matFiles = dir(fullfile(input_folder,'*.mat'));
Tmin = 2; %seconds
Tmax = 5; %seconds

for fileID = 1:length(matFiles)
    rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
    if (isfield(rawData,'segment'))
        rawData = rawData.segment;
    elseif (isfield(rawData, 'rawData'))
        rawData = rawData.rawData;
    end
    if (isfield(rawData, "c02_left"))
        rawData.vehicleOppositeSideDisturbing = double(abs(rawData.c02_left) < 0.5);
        counter = 0;
        startIdx = 0;
        dt = mean(diff(rawData.q_T0));
        for i=1:length(rawData.vehicleOppositeSideDisturbing)
            if (rawData.vehicleOppositeSideDisturbing(i) == 1 && counter == 0)
                startIdx = i;
                counter = counter + dt;
            elseif (rawData.vehicleOppositeSideDisturbing(i) == 1)
                counter = counter + dt;
            else
                if ((counter < Tmin || counter > Tmax) && startIdx > 0)
                    rawData.vehicleOppositeSideDisturbing(startIdx:i) = 0;
                elseif (startIdx > 0)
                    rawData.vehicleOppositeSideDisturbing(startIdx:i) = counter;
                end
                counter = 0;
                startIdx = 0;
            end
        end
        save(fullfile(output_folder,strcat(matFiles(fileID).name(1:end-4),'_withTrafficInfo.mat')),'rawData');
    else
        disp('Missing opposite lane info in the following file:');
        disp(matFiles(fileID).name);
    end
end