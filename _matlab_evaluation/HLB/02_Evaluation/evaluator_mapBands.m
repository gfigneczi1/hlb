function evaluator_mapBands(segment,config)
% This evaluator is responsible for generating the reference_polygon mat
% file (under Inputs folder) which contains polygons, which describe the
% valid areas around the globe. An area is valid, if a road which falls
% into it is useful for evaluating.
% Example: road 62/Hungary between Dunaujvaros and Szekesfehervar. The
% polygon approximately falls to the lane edges of this road. The polygon
% exludes roundabouts or side roads, or sections which are useless in
% certain aspect.

% Segmentor: segmentor_mapBand

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-08

dLat = 2e-4;
dLong = 2e-4;
%% initialize start/stop polygon points for each road segments
% 31 road - reference S-curve
roadPoint = [47.44727500358191, 19.382351244924102];

Longv_start_31_ref = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_start_31_ref = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_start_31_ref = Longv_start_31_ref * 40075000 .* cos(Latv_start_31_ref*pi()/180) / 360;
yv_start_31_ref = Latv_start_31_ref * 111.32 * 1000;
roadPoint = [47.44194776636965, 19.400457497390107]; 
Longv_stop_31_ref = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_stop_31_ref = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_stop_31_ref = Longv_stop_31_ref * 40075000 .* cos(Latv_stop_31_ref*pi()/180) / 360;
yv_stop_31_ref = Latv_stop_31_ref * 111.32 * 1000;

% M31 - M0 to Maglód
roadPoint = [47.45958159037517, 19.343089470216757];

% always: roadPoint(1) - dLat / roadPoint(1) - dLat / roadPoint(1) + dLat /
% roadPoints(1) + dLat (starting from box bottom left, going clockwise
% then for LONG: roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong,
% roadPoint(2) - dLong

Longv_start_31_1 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_start_31_1 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_start_31_1 = Longv_start_31_1 * 40075000 .* cos(Latv_start_31_1*pi()/180) / 360;
yv_start_31_1 = Latv_start_31_1 * 111.32 * 1000;
roadPoint = [47.45146490583005, 19.358584360288816]; 
Longv_stop_31_1 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_stop_31_1 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_stop_31_1 = Longv_stop_31_1 * 40075000 .* cos(Latv_stop_31_1*pi()/180) / 360;
yv_stop_31_1 = Latv_stop_31_1 * 111.32 * 1000;

% 31 - Maglód to Mende
roadPoint = [47.44777134561812, 19.376820859974515];
Longv_start_31_2 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_start_31_2 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_start_31_2 = Longv_start_31_2 * 40075000 .* cos(Latv_start_31_2*pi()/180) / 360;
yv_start_31_2 = Latv_start_31_2 * 111.32 * 1000;

roadPoint = [47.42647856880319, 19.437365872376215];
Longv_stop_31_2 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_stop_31_2 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_stop_31_2 = Longv_stop_31_2 * 40075000 .* cos(Latv_stop_31_2*pi()/180) / 360;
yv_stop_31_2 = Latv_stop_31_2 * 111.32 * 1000;

% 31 - Mende - Sülysáp
roadPoint = [47.43596013773291, 19.466462505285392];
Longv_start_31_3 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_start_31_3 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_start_31_3 = Longv_start_31_3 * 40075000 .* cos(Latv_start_31_3*pi()/180) / 360;
yv_start_31_3 = Latv_start_31_3 * 111.32 * 1000;

roadPoint = [47.44689897838279, 19.521105062050783];
Longv_stop_31_3 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_stop_31_3 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_stop_31_3 = Longv_stop_31_3 * 40075000 .* cos(Latv_stop_31_3*pi()/180) / 360;
yv_stop_31_3 = Latv_stop_31_3 * 111.32 * 1000;

% 31 - Tapioszecso - Szentmartonkata
roadPoint = [47.458336558176086, 19.60221522752861];
Longv_start_31_4 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_start_31_4 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_start_31_4 = Longv_start_31_4 * 40075000 .* cos(Latv_start_31_4*pi()/180) / 360;
yv_start_31_4 = Latv_start_31_4 * 111.32 * 1000;

roadPoint = [47.46294366335314, 19.65318283693067];
Longv_stop_31_4 = [roadPoint(2) - dLong, roadPoint(2) + dLong, roadPoint(2) + dLong , roadPoint(2) - dLong ];
Latv_stop_31_4 = [roadPoint(1) - dLat, roadPoint(1) - dLat, roadPoint(1) + dLat, roadPoint(1) + dLat];
xv_stop_31_4 = Longv_stop_31_4 * 40075000 .* cos(Latv_stop_31_4*pi()/180) / 360;
yv_stop_31_4 = Latv_stop_31_4 * 111.32 * 1000;

% Zala highway
Longv_start_hw = [16.842256, 16.842370800476, 16.842385 ,16.842286 ];
Latv_start_hw = [46.897053,46.89707452274891 , 46.896985 ,46.896969 ];
xv_start_hw = Longv_start_hw * 40075000 .* cos(Latv_start_hw*pi()/180) / 360;
yv_start_hw = Latv_start_hw * 111.32 * 1000;
Longv_stop_hw = [16.83938259524578, 16.839409417333645, 16.839520728998302, 16.83949055414945];
Latv_stop_hw = [46.892269684315906, 46.892073558600906, 46.89204606428449, 46.892278849051394];
xv_stop_hw = Longv_stop_hw * 40075000 .* cos(Latv_stop_hw*pi()/180) / 360;
yv_stop_hw = Latv_stop_hw * 111.32 * 1000;

% Zala Rural road southern part
Longv_start_Rural_south = [16.846885, 16.846885, 16.846895, 16.846895];
Latv_start_Rural_south = [46.88815, 46.88818, 46.88815, 46.88818];
xv_start_Rural_south = Longv_start_Rural_south * 40075000 .* cos(Latv_start_Rural_south*pi()/180) / 360;
yv_start_Rural_south = Latv_start_Rural_south * 111.32 * 1000;
Longv_stop_Rural_south = [16.84123, 16.84143, 16.841083, 16.841338];
Latv_stop_Rural_south = [46.88467, 46.88508, 46.884401, 46.884486];
xv_stop_Rural_south = Longv_stop_Rural_south * 40075000 .* cos(Latv_stop_Rural_south*pi()/180) / 360;
yv_stop_Rural_south = Latv_stop_Rural_south * 111.32 * 1000;

 % Zala second after roundabout
Longv_start_Rural_south = [16.83936, 16.839242, 16.839518, 16.839643];
Latv_start_Rural_south = [46.88431, 46.884209, 46.88408, 46.884188];
xv_start_Rural_south_2 = Longv_start_Rural_south * 40075000 .* cos(Latv_start_Rural_south*pi()/180) / 360;
yv_start_Rural_south_2 = Latv_start_Rural_south * 111.32 * 1000;
Longv_stop_Rural_south = [16.846924, 16.846926, 16.847302, 16.847297];
Latv_stop_Rural_south = [46.88596, 46.88586, 46.885851, 46.885928];
xv_stop_Rural_south_2 = Longv_stop_Rural_south * 40075000 .* cos(Latv_stop_Rural_south*pi()/180) / 360;
yv_stop_Rural_south_2 = Latv_stop_Rural_south * 111.32 * 1000;
% M6 highway
Longv_start_M6 = [18.969015, 18.969585, 18.970099, 18.969418];
Lateralv_start_M6 = [47.394933, 47.395006, 47.393403, 47.393359];
xv_start_M6 = Longv_start_M6 * 40075000 .* cos(Lateralv_start_M6*pi()/180) / 360;
yv_start_M6 = Lateralv_start_M6 * 111.32 * 1000;
Longv_end_M6 = [18.885888, 18.886410, 18.887478, 18.886905];
Lateralv_end_M6 = [46.997369, 46.997335, 47.000746, 47.000817];
xv_end_M6 = Longv_end_M6 * 40075000 .* cos(Lateralv_end_M6*pi()/180) / 360;
yv_end_M6 = Lateralv_end_M6 * 111.32 * 1000;
% M7 highway
Longv_start_M7 = [18.485658, 18.490669, 18.490156, 18.485244];
Lateralv_start_M7 = [47.180092, 47.182135, 47.182611, 47.180513];
xv_start_M7 = Longv_start_M7 * 40075000 .* cos(Lateralv_start_M7*pi()/180) / 360;
yv_start_M7 = Lateralv_start_M7 * 111.32 * 1000;
Longv_end_M7 = [18.885590, 18.885171, 18.882838, 18.883466];
Lateralv_end_M7 = [47.421519, 47.421641, 47.418708, 47.418563];
xv_end_M7 = Longv_end_M7 * 40075000 .* cos(Lateralv_end_M7*pi()/180) / 360;
yv_end_M7 = Lateralv_end_M7 * 111.32 * 1000;
% 62 ROAD (1,2,3 roundabouts)
% 62 start to first roundabout start
Longv_start_62_1 = [18.877465, 18.877270, 18.872280, 18.872656];
Lateralv_start_62_1 = [46.999212, 46.999134, 47.002367, 47.002565];
xv_start_62_1 = Longv_start_62_1 * 40075000 .* cos(Lateralv_start_62_1*pi()/180) / 360;
yv_start_62_1 = Lateralv_start_62_1 * 111.32 * 1000;
Longv_end_62_1 = [18.809303, 18.809013, 18.810314, 18.810873];
Lateralv_end_62_1 = [47.056201, 47.056060, 47.054717, 47.054912];
xv_end_62_1 = Longv_end_62_1 * 40075000 .* cos(Lateralv_end_62_1*pi()/180) / 360;
yv_end_62_1 = Lateralv_end_62_1 * 111.32 * 1000;
% 62 first roundabout end to second roundabout start
Longv_start_62_2 = [18.808534, 18.808129, 18.806381, 18.806879];
Lateralv_start_62_2 = [47.056913, 47.056772, 47.058180, 47.058416];
xv_start_62_2 = Longv_start_62_2 * 40075000 .* cos(Lateralv_start_62_2*pi()/180) / 360;
yv_start_62_2 = Lateralv_start_62_2 * 111.32 * 1000;
Longv_end_62_2 = [18.629050, 18.629088, 18.631389, 18.631362];
Lateralv_end_62_2 = [47.105438, 47.105197, 47.105271, 47.105607];
xv_end_62_2 = Longv_end_62_2 * 40075000 .* cos(Lateralv_end_62_2*pi()/180) / 360;
yv_end_62_2 = Lateralv_end_62_2 * 111.32 * 1000;
% 62 second roundabout end to third roundabout start
Longv_start_62_3 = [18.627498, 18.627538, 18.624864, 18.624675];
Lateralv_start_62_3 = [47.105330, 47.105049, 47.104552, 47.104910];
xv_start_62_3 = Longv_start_62_3 * 40075000 .* cos(Lateralv_start_62_3*pi()/180) / 360;
yv_start_62_3 = Lateralv_start_62_3 * 111.32 * 1000;
Longv_end_62_3 = [18.594219, 18.594444, 18.596468, 18.596312];
Lateralv_end_62_3 = [47.094848, 47.094585, 47.095236, 47.095484];
xv_end_62_3 = Longv_end_62_3 * 40075000 .* cos(Lateralv_end_62_3*pi()/180) / 360;
yv_end_62_3 = Lateralv_end_62_3 * 111.32 * 1000;
% 62 third roundabout end to 62 end
Longv_start_62_4 = [18.591771, 18.591894, 18.588957, 18.588912];
Lateralv_start_62_4 = [47.094234, 47.093898, 47.093502, 47.093805];
xv_start_62_4 = Longv_start_62_4 * 40075000 .* cos(Lateralv_start_62_4*pi()/180) / 360;
yv_start_62_4 = Lateralv_start_62_4 * 111.32 * 1000;
Longv_end_62_4 = [18.485612, 18.485365, 18.488010, 18.488308];
Lateralv_end_62_4 = [47.175769, 47.175519, 47.174316, 47.174566];
xv_end_62_4 = Longv_end_62_4 * 40075000 .* cos(Lateralv_end_62_4*pi()/180) / 360;
yv_end_62_4 = Lateralv_end_62_4 * 111.32 * 1000;
% M4 from Budapest to Szolnok
Longv_start_M4_1 = [19.32422530588145, 19.328299024032418, 19.328649120980092, 19.32417708060526];
Latv_start_M4_1 = [47.404097382870084, 47.403809353951075,47.40402635448332 ,47.40392738262679];
Longv_end_M4_1 = [19.808199164260778, 19.809629777015214, 19.809267563585585, 19.807405080338537];
Latv_end_M4_1 = [47.19700374314316, 47.19675851120928,47.197029417179 ,47.19735387795167];
xv_start_M4 = Longv_start_M4_1 * 40075000 .* cos(Latv_start_M4_1*pi()/180) / 360;
yv_start_M4 = Latv_start_M4_1 * 111.32 * 1000;
xv_end_M4 = Longv_end_M4_1 * 40075000 .* cos(Latv_end_M4_1*pi()/180) / 360;
yv_end_M4 = Latv_end_M4_1 * 111.32 * 1000;
% initializing points to check
X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
Y_abs = segment.LatPos_abs * 111.32*1000;
%% start checking which point is within the polygon
% fields should start with the road type (M for highway, r for country road)
road_segments.M6.in_start_M6 = inpolygon(X_abs, Y_abs, xv_start_M6, yv_start_M6);
road_segments.M6.in_end_M6 = inpolygon(X_abs, Y_abs, xv_end_M6, yv_end_M6);

road_segments.r62_1.in_start_62_1 = inpolygon(X_abs, Y_abs, xv_start_62_1, yv_start_62_1);
road_segments.r62_1.in_end_62_1 = inpolygon(X_abs, Y_abs, xv_end_62_1, yv_end_62_1);
road_segments.r62_2.in_start_62_2 = inpolygon(X_abs, Y_abs, xv_start_62_2, yv_start_62_2);
road_segments.r62_2.in_end_62_2 = inpolygon(X_abs, Y_abs, xv_end_62_2, yv_end_62_2);
road_segments.r62_3.in_start_62_3 = inpolygon(X_abs, Y_abs, xv_start_62_3, yv_start_62_3);
road_segments.r62_3.in_end_62_3 = inpolygon(X_abs, Y_abs, xv_end_62_3, yv_end_62_3);
road_segments.r62_4.in_start_62_4 = inpolygon(X_abs, Y_abs, xv_start_62_4, yv_start_62_4);
road_segments.r62_4.in_end_62_4 = inpolygon(X_abs, Y_abs, xv_end_62_4, yv_end_62_4);

road_segments.r31_1.in_start_31_1 = inpolygon(X_abs, Y_abs, xv_start_31_1, yv_start_31_1);
road_segments.r31_1.in_end_31_1 = inpolygon(X_abs, Y_abs, xv_stop_31_1, yv_stop_31_1);
road_segments.r31_2.in_start_31_2 = inpolygon(X_abs, Y_abs, xv_start_31_2, yv_start_31_2);
road_segments.r31_2.in_end_31_2 = inpolygon(X_abs, Y_abs, xv_stop_31_2, yv_stop_31_2);
road_segments.r31_3.in_start_31_3 = inpolygon(X_abs, Y_abs, xv_start_31_3, yv_start_31_3);
road_segments.r31_3.in_end_31_3 = inpolygon(X_abs, Y_abs, xv_stop_31_3, yv_stop_31_3);
road_segments.r31_4.in_start_31_4 = inpolygon(X_abs, Y_abs, xv_start_31_4, yv_start_31_4);
road_segments.r31_4.in_end_31_4 = inpolygon(X_abs, Y_abs, xv_stop_31_4, yv_stop_31_4);

road_segments.ref31.in_start_31 = inpolygon(X_abs, Y_abs, xv_start_31_ref, yv_start_31_ref);
road_segments.ref31.in_end_31 = inpolygon(X_abs, Y_abs, xv_stop_31_ref, yv_stop_31_ref);

road_segments.M7.in_start_M7 = inpolygon(X_abs, Y_abs, xv_start_M7, yv_start_M7);
road_segments.M7.in_end_M7 = inpolygon(X_abs, Y_abs, xv_end_M7, yv_end_M7);

road_segments.M4.in_start_M4 = inpolygon(X_abs, Y_abs, xv_start_M4, yv_start_M4);
road_segments.M4.in_end_M4 = inpolygon(X_abs, Y_abs, xv_end_M4, yv_end_M4);

road_segments.ZZ_hw.in_start_hw = inpolygon(X_abs, Y_abs, xv_start_hw, yv_start_hw);
road_segments.ZZ_hw.in_end_hw = inpolygon(X_abs, Y_abs, xv_stop_hw, yv_stop_hw);
road_segments.ZZ_rural_south.in_start_rural_south = inpolygon(X_abs, Y_abs, xv_start_Rural_south, yv_start_Rural_south);
road_segments.ZZ_rural_south.in_end_rural_south = inpolygon(X_abs, Y_abs, xv_stop_Rural_south, yv_stop_Rural_south);
road_segments.ZZ_rural_south_2.in_start_rural_south_2 = inpolygon(X_abs, Y_abs, xv_start_Rural_south_2, yv_start_Rural_south_2);
road_segments.ZZ_rural_south_2.in_end_rural_south_2 = inpolygon(X_abs, Y_abs, xv_stop_Rural_south_2, yv_stop_Rural_south_2);
%% creating polygon for road segments with bandwidth
fn_rs = fieldnames(road_segments);
bandwidth = 25; %[m]
dy_right = -bandwidth * 1/4; %[m]
dy_left = bandwidth * 3/4; %[m]
for i=1:numel(fn_rs)
    road_segment = road_segments.(fn_rs{i});
    fn_s = fieldnames(road_segment);
    start_poly = fn_s{1};
    end_poly = fn_s{2};
    if numel(X_abs(road_segment.(start_poly))) > 0
        points_in_start = X_abs(road_segment.(start_poly));
        points_in_end = X_abs(road_segment.(end_poly));
        % if the measurement started later then the first checkpoint or
        % ended sooner then the last checkpoint, use the first/last
        % coordinates
        if isempty(points_in_start) && start_poly == "in_start_M6"
            points_in_start = X_abs(1);
        elseif isempty(points_in_end) && end_poly == "in_end_M7"
            points_in_end = X_abs(end);
        end
        % index generation for slicing
        if (isempty(points_in_start))
            slice_start = 1;
        else
            slice_start = find(X_abs == points_in_start(1));
        end
        if (isempty(points_in_end))
            slice_end = length(X_abs);
        else
            slice_end = find(X_abs == points_in_end(end));
        end
        if (slice_end(1) < slice_start(1))
            temp = slice_start(1);
            slice_start(1) = slice_end(1);
            slice_end(1) = temp;
            clear temp;
        end
        % calculating offset polygon points
        theta = segment.theta_calc(slice_start(1):slice_end(end));
        X_poly_right = -sin(theta)*dy_right + X_abs(slice_start(1):slice_end(end));
        Y_poly_right = cos(theta)*dy_right + Y_abs(slice_start(1):slice_end(end));
        X_poly_left = -sin(theta)*dy_left + X_abs(slice_start(1):slice_end(end));
        Y_poly_left = cos(theta)*dy_left + Y_abs(slice_start(1):slice_end(end));
        X_poly = [X_poly_right; X_poly_left];
        Y_poly = [Y_poly_right; Y_poly_left];
        % reducing the number of points then save it by X,Y coordinates
        % 0.00000098 is a great parameter for highways, but little bit
        % small for country roads or test roads
        reduced_polygon = reducepoly([X_poly, Y_poly], 0.00000098);
        reference_polygon.(fn_rs{i}).X_poly = reduced_polygon(:,1);
        reference_polygon.(fn_rs{i}).Y_poly = reduced_polygon(:,2);
    end
end
%% reading the existing reference polygon first
reference_polygon_in = load(fullfile('Inputs','reference_polygon.mat'));
%% saving the created polygon to a mat file
fn_rs = fieldnames(reference_polygon);
for i=1:length(fn_rs)
    reference_polygon_in.(fn_rs{i}) = reference_polygon.(fn_rs{i});
end
path = fullfile('Inputs','reference_polygon.mat');
save(path, '-struct', 'reference_polygon_in');
end

