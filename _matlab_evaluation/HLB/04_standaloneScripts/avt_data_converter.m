clear; close all; clc;

path = "\\abtvdfs2.de.bosch.com\ismdfs\ilm\abt\AVT_Datasets\CAN_Exports\Cadillac_CT6";
mappingFile = "\\abtvdfs2.de.bosch.com\ismdfs\ilm\abt\AVT_Datasets\trip_keys\2023-12-06_Archive_check_update.xlsx";
outputFolder = "C:\database\AVT\Cadillac";

necessarySignals = ["LeftLaneMrkLatPos_0x346_0", "RightLaneMrkLatPos_0x347_0", "LaneSenseModeState_0x345_0", ...
    "ts_s_0", ...
    "PPSLong_0x262_100", "PPSLat_0x261_100", "PPSHeading_0x260_100", ...
    "LongitudinalVelocity_0x344_20", "LateralVelocity_0x344_20", "VehicleSpeedAvgDrvn_0x3E9_100", ...
    "InertialLongitudinalAcc_0x340_10", "InertialLateralAcc_0x340_10", ...
    "LonAccSnsVal_0x1FC_50", "LatAccSnsVal_0x1FC_50", ...
    "IMUProtYawRt_0x34C_50", "YawRate_0x1E9_20", ...
    "ACCTorqueCmdAxleTorqueRequest_0x2CB_0", ...
    "FrwayRoadTypInfo_0x150_100", "LKASteeringCmdActive_0x152_10", "LKASMode_0x152_10", ...
    "ACCactive_0x2CB_0", "BrakePedalPos_0x0BE_12", "AcceleratorPedal_0x1A1_25", ...
    "LKADriverAppldTrq_0x164_10", "LKATBDTorque_0x164_10", ...
    "SteeringWheelAngle_0x1E5_10", "SteeringWheelAngleGradient_0x1E5_10"];

% data details and units, resolution etc.
% LaneMrkLatPos: meters
% PPSLong/Lat: degrees
% LongitudinalVelocity_0x344_20 (and lat): m/s
% InertialLongitudinalAcc_0x340_10 (and lat): m/s^2
% LonAccSnsVal_0x1FC_50 (and lat): g
% IMUProtYawRt_0x34C_50: deg/s
% YawRate_0x1E9_20: deg/s
% ACCTorqueCmdAxleTorqueRequest_0x2CB_0: Nm
% FrwayRoadTypInfo_0x150_100: 0 = NonFreeway, 1 = Freeway
% LKASMode_0x152_10: 0 none, 1< some kind of control
% BrakePedalPos_0x0BE_12: %
% AcceleratorPedal_0x1A1_25: %
% LKADriverAppldTrq_0x164_10: Nm
% LKATBDTorque_0x164_10: Nm
% SteeringWheelAngle_0x1E5_10: deg
% SteeringWheelAngleGradient_0x1E5_10: deg/sec    

map = readtable(mappingFile);
tripNames = map.trip_label;

csvFiles = dir(fullfile(path, "*.csv"));

for i=1:length(csvFiles)
    data = readtable(fullfile(csvFiles(i).folder, csvFiles(i).name));
    preFilterCriteriaLaneOffset = ~all(data.LeftLaneMrkLatPos_0x346_0==0);
    preFilterCriteriaAcceleration = all(data.ACCactive_0x2CB_0==0);
    a = find(csvFiles(i).name=='_',3);    
    dr = map.subject_nr(find(convertCharsToStrings(tripNames)==convertCharsToStrings(csvFiles(i).name(1:a(3)-1))));
    tripId = map.trip_id(find(convertCharsToStrings(tripNames)==convertCharsToStrings(csvFiles(i).name(1:a(3)-1))));
    year = map.year(find(convertCharsToStrings(tripNames)==convertCharsToStrings(csvFiles(i).name(1:a(3)-1))));
    month = map.month(find(convertCharsToStrings(tripNames)==convertCharsToStrings(csvFiles(i).name(1:a(3)-1))));
    day = map.day(find(convertCharsToStrings(tripNames)==convertCharsToStrings(csvFiles(i).name(1:a(3)-1))));
    
    if (preFilterCriteriaAcceleration)    
       for j=1:numel(necessarySignals)
           dataMasked.(necessarySignals(j)) = data.(necessarySignals(j));
       end
       save(fullfile(outputFolder, strcat(num2str(dr),'_', num2str(tripId), '_', num2str(year), num2str(month), num2str(day), '.mat')), '-struct', 'dataMasked');
    end
    
    clc;
    fprintf("Files processed:\n");
    fprintf(strcat(num2str(i), "/", num2str(length(csvFiles))));
end