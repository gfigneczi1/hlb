#################################################################################################
#                                                                                               #
#   ATTENTION: This file has to be run with a Python environment of Python 3.10 or higher       #
#   This is because it is dependent on geojson and shapely                                      #
#                                                                                               #
#################################################################################################
import geopandas as gpd
from shapely.geometry import LineString, Point
import pandas as pd
import numpy as np
import pathlib
import os
import sys
import json
import math
import time

LONGITUDE_BOUND_IN_METERS = 1000
LATITUDE_BOUND_IN_METERS = 1000

CURVATURE_TIME_INTERVAL = 60

def convert_position_to_meters(latitude, longitude):
    '''returns longitude and latitude in meters instead of degrees

    keyword arguments:
    latitude -- latitude of the point we want to convert 
    longitude -- longitude of the point we want to convert
    '''
    meter_latitude = 111320 * latitude
    meter_longitude = longitude * 40075000 * math.cos(math.radians(latitude))/360
    return meter_longitude, meter_latitude

def find_closest_point(latitude, longitude, coords):
    '''returns the closest point to a coordinate on a given street line

    keyword arguments:
    latitude -- latitude of the point we want to find the closest point on the given street to
    longitude -- longitude of the point we want to find the closest point on the given street to
    coords -- list of coordinates of all the points on the street in question
    '''
    meter_longitude, meter_latitude = convert_position_to_meters(latitude, longitude)
    smallest_distance = 10000000000
    smallest_index = -1
    for i, coo in enumerate(coords):
        meter_road_longitude, meter_road_latitude = convert_position_to_meters(coo[0], coo[1])
        a = abs(meter_longitude - meter_road_longitude)
        b = abs(meter_latitude - meter_road_longitude)
        c = math.sqrt(math.pow(a, 2) + math.pow(b, 2))
        if c < smallest_distance:
            smallest_distance = c
            smallest_index = i
    return smallest_index

def extract_guidance_info(guidance_dict, attribute, attribute_detail):
    '''    returns the information about the attribute including the coordinates in a dictionary
    special case: returns 'NaN' if the dictionary has no info on the guidance attribute

    keyword arguments:
    guidance_dict -- dictionary extracted from a Street Line element with the 'guidance' attribute
    attribute -- The attribute that is extracted (TRAFFIC_LIGHTS, PEDESTRIAN_CROSSINGS, STOP_SIGNS)
    attribute_detail -- The attribute that is extracted more specifically
    '''
    generic_dict = {}
    if not isinstance(guidance_dict, type(generic_dict)):
        return str(math.nan)
    elif attribute in guidance_dict.keys():
        sub_dict = guidance_dict[attribute]
        if not isinstance(sub_dict, type(generic_dict)):
            return math.nan
        elif attribute_detail in sub_dict.keys():
            if attribute == 'TRAFFIC_LIGHTS':
                print(guidance_dict[attribute][attribute_detail])
            return guidance_dict[attribute][attribute_detail]
        else:
            return guidance_dict[attribute]
    else:
        return math.nan

def extract_guidance_overall_info(guidance_dict, attribute):
    '''returns the information about the attribute including the coordinates in a dictionary
    special case: returns 'NaN' if the dictionary has no info on the guidance attribute

    keyword arguments:
    guidance_dict -- dictionary extracted from a Street Line element with the 'guidance' attribute
    attribute -- The attribute that is extracted (TRAFFIC_LIGHTS, PEDESTRIAN_CROSSINGS, STOP_SIGNS)
    '''
    generic_dict = {}
    if not isinstance(guidance_dict, type(generic_dict)):
        return math.nan
    elif attribute in guidance_dict.keys():
        sub_dict = guidance_dict[attribute]
        return sub_dict
    else:
        return math.nan

def extract_from_link_type(gdf, nearest_index, search_link_type):
    '''returns True if the link type corresponds to the link type in question

    keyword arguments:
    gdf -- the GeoPandas GeoDataFrame with the extracted geographic information about the trip
    nearest_index -- the index of the closest road element
    search_link_type -- the type of link in question (options: ramp, roundabout)
    '''
    link_type = gdf['linkType'][nearest_index]
    if link_type == search_link_type:
        return True
    else:
        return False

def within_bounds(meter_longitude, meter_latitude, road_element_nearby_bounds_meters):
    '''returns True if the coordinates of a point are within the bounds of a box around a road element like a traffic light or a pedestrian crossing

    keyword arguments:
    meter_longitude -- longitude of the gps coordinate point converted to meters
    meter_latitude -- latitude of the gps coordinate point converted to meters
    road_element_nearby_bounds_meters -- list where the elements are the coordinates in meters of the box around a road element,
        an element is formatted as follows: [west_bound, east_bound, south_bound, north_bound]
    '''
    for road_element in road_element_nearby_bounds_meters:
        if meter_longitude > road_element[0] and meter_longitude < road_element[1]:
            if meter_latitude > road_element[2] and meter_latitude < road_element[3]:
                return True
    return False

def extract_bound_coordinates(road_element_list):
    '''returns list where the elements are the coordinates in meters of the box around a road element
    The list is formatted as follows [west_bound, east_bound, south_bound, north_bound]

    keyword arguments:
    road_element_list -- list of dictionaries with the guidance attribute information (traffic_light, pedestrian_crossing)
    '''
    road_element_nearby_bounds = []
    for i, element in enumerate(road_element_list):
        if isinstance(element, dict):
            type_string = str(type(element))
            if 'dict' in type_string:
                validity_dict = element['_validity']
            elif 'list' in type_string:
                validity_dict = element[0]['_validity']
            else:
                print(type(element))
            meter_longitude, meter_latitude = convert_position_to_meters(longitude=validity_dict['coordinates'][0][0], latitude=validity_dict['coordinates'][0][1])
            road_element_nearby_bounds.append([meter_longitude-LONGITUDE_BOUND_IN_METERS, meter_longitude+LONGITUDE_BOUND_IN_METERS, meter_latitude-LATITUDE_BOUND_IN_METERS, meter_latitude+LATITUDE_BOUND_IN_METERS])
    return road_element_nearby_bounds

def extract_stop_sign_bound_coordinates(warning_signs, warning_signs_meanings):
    '''returns list where the elements are the coordinates in meters of the box around a stop sign
    the list is formatted as follows: [west_bound, east_bound, south_bound, north_bound]

    keyword arguments:
    warning_signs -- list of dictionaries including warning sign information (e.g. coordinates)
    warning_sign_meanings -- list of the meanings of the warning signs (e.g. stop sign)
    '''
    stop_sign_bounds_meters = []
    for i, (meaning, sign) in enumerate(zip(warning_signs_meanings, warning_signs)):
        if isinstance(sign, dict):
            type_string = str(type(sign))
            if 'dict' in type_string:
                validity_dict = sign['_validity']
            elif 'list' in type_string:
                validity_dict = sign[0]['_validity']
            else:
                print(type(element))
        if 'STOP' in str(meaning):
            meter_longitude, meter_latitude = convert_position_to_meters(validity_dict['coordinates'][0][0], validity_dict['coordinates'][0][1])
            stop_sign_bounds_meters.append([meter_longitude-LONGITUDE_BOUND_IN_METERS, meter_longitude+LONGITUDE_BOUND_IN_METERS, meter_latitude-LATITUDE_BOUND_IN_METERS, meter_latitude+LATITUDE_BOUND_IN_METERS])
    return stop_sign_bounds_meters

def convert_slope_value_to_percentage(slope_value):
    '''returns the value of the slope in %
    special case: returns 'NaN' if the slope is an invalid value, not accounted for in the NDS encoding

    keyword arguments:
    slope_value -- the value is encoded in the NDS encoding (https://doc.nds-association.org/version/2.5.4/code/nds.common.flexattr.attrdefs.Slope)
    '''
    if slope_value == 126:
        return 30.1
    elif slope_value == 125:
        return 27.5
    elif slope_value > -124 and slope_value < 124:
        return slope_value * 0.2
    else:
        print("ATTENTION: invalid slope value")
        return math.nan

def get_slope_value(street_index, point_index, gdf):
    '''returns returns the slope value in percentage after extracting it 
    special case: returns 'NaN' if it cannot find the slope

    keyword arguments:
    street_index -- index of the Line element representing the street in question in the GeoPandas GeoDataFrame        
    point_index -- index of the point inside the LineString representing the street in question (every LineString is made out of points)
    gdf -- the GeoPandas GeoDataFrame with the extracted geographic information about the trip
    '''
    if 'dict' in str(type(gdf["adas"][street_index])):
        slope_point_list = gdf["adas"][street_index]['SLOPE_ARRAY']['slopeArray']['slopePoint']
        for i in range(len(slope_point_list)):
            slope_point_index = slope_point_list[i]['point']['numVx4']
            if slope_point_index == point_index:
                slope_value = slope_point_list[i]['slope']
                slope_value_percentage = convert_slope_value_to_percentage(slope_value)
                return slope_value_percentage
    return math.nan

def get_curvature_value(point1_latitude, point1_longitude, point2_latitude, point2_longitude, point3_latitude, point3_longitude):
    '''returns a curvature calculated from gps points, the time distance between the sampled gps points
    This calculation is based on the 3 point triangle formula to calculate the curve radius
    
    keyword arguments:
    point1_latitude -- latitude of the first point on the curve
    point1_longitude -- longitude of the first point on the curve
    point2_latitude -- latitude of the middle point on the curve
    point2_longitude -- longitude of the middle point on the curve
    point3_latitude -- latitude of the last point on the curve
    point3_longitude -- longitude of the last point on  the curve
    '''
    #the point latitudes and longitudes are converted to meters
    point1_meter_longitude, point1_meter_latitude = convert_position_to_meters(point1_latitude, point1_longitude)
    point2_meter_longitude, point2_meter_latitude = convert_position_to_meters(point2_latitude, point2_longitude)
    point3_meter_longitude, point3_meter_latitude = convert_position_to_meters(point3_latitude, point3_longitude)
    # distances between the points in meters are calculated here
    a = math.sqrt((point1_meter_latitude - point2_meter_latitude)*(point1_meter_latitude - point2_meter_latitude) + (point1_longitude - point2_longitude)*(point1_longitude - point2_longitude))
    b = math.sqrt((point1_meter_latitude - point3_meter_latitude)*(point1_meter_latitude - point3_meter_latitude) + (point1_longitude - point3_longitude)*(point1_longitude - point3_longitude))
    c = math.sqrt((point3_meter_latitude - point2_meter_latitude)*(point3_meter_latitude - point2_meter_latitude) + (point3_longitude - point2_longitude)*(point3_longitude - point2_longitude))
    s = (a+b+c)/2
    triangle_area_between_points = math.sqrt(abs(s*(s-a)*(s-b)*(s-c))) #area of the triangle between the three points
    if a == 0 or c == 0 or b == 0:
        return math.nan
    else:
        C = 4*triangle_area_between_points/(a*b*c)
        return C

def get_map_element_instances(gdf, latitudes, longitudes):
    '''returns the instances of the various map elements extracted from the NDS map as time sequences
    Outputted time sequences
        - NDS location information: nearest_indexes, nearest_elements
        - binary status columns: motorway, urban, ramps, roundabouts
        - non-binary states: link_types, warning_signs_meanings
        - dictionaries used later: overall_traffic_lights, pedestrian_crossings, warning_signs
        - int values:  averageSpeeds, speed_limits, advisory_speed_limits
        - float values: slope_values_in_percentages, curvatures
    keyword arguments:
    gdf -- the GeoPandas GeoDataFrame with the extracted geographic information about the trip
    latitudes -- list of the latitude time sequence points found in the recording
    longitudes -- list of the longitude time sequence points found in the recording
    '''
    nearest_indexes = []
    nearest_elements = []
    #broad road types
    urban = []
    motorway = []
    link_types = []
    averageSpeeds = []
    #guidance info
    speed_limits = []
    advisory_speed_limits = []
    traffic_lights = []
    pedestrian_crossings = []
    warning_signs = []
    warning_signs_meanings = []
    roundabouts = []
    ramps = []
    overall_traffic_lights = []
    #slope
    slope_values_in_percentages = []
    #curvature
    curvatures = []
    for i, (latitude, longitude) in enumerate(zip(latitudes, longitudes)):
        gps_point = Point(longitude, latitude)
        #get nearest geometrical element to position
        nearest_geometry = gdf.sindex.nearest(gps_point, return_all=False)
        nearest_index = int(nearest_geometry[1])
        nearest_indexes.append(nearest_index)
        nearest_element = gdf["id"][nearest_index]
        nearest_elements.append(nearest_element)
        #to get slope
        nearest_element_coordinates = list(gdf["geometry"][nearest_index].coords)
        closest_point_index = find_closest_point(latitude=latitude, longitude=longitude, coords=nearest_element_coordinates)
        nearest_road_point = nearest_element_coordinates[closest_point_index]
        slope_values_in_percentages.append(get_slope_value(nearest_index, closest_point_index, gdf))
        #to calculate curve
        if i > CURVATURE_TIME_INTERVAL and i < len(latitudes)-CURVATURE_TIME_INTERVAL:
            curve = get_curvature_value(latitudes[i-CURVATURE_TIME_INTERVAL], longitudes[i-CURVATURE_TIME_INTERVAL], latitude, longitude, latitudes[i+CURVATURE_TIME_INTERVAL], longitudes[i+CURVATURE_TIME_INTERVAL])
        else:
            curve = math.nan
        curvatures.append(curve)
        #general map infos
        urban.append(gdf['urban'][nearest_index])
        motorway.append(gdf['motorway'][nearest_index])
        link_types.append(gdf['linkType'][nearest_index])
        roundabouts.append(extract_from_link_type(gdf, nearest_index, 'ROUNDABOUT'))
        ramps.append(extract_from_link_type(gdf, nearest_index, 'RAMP'))
        averageSpeeds.append(gdf['averageSpeed'][nearest_index])
        #these attributes are hidden in the guidance attribute
        guidance_dict = gdf["guidance"][nearest_index]
        speed_limits.append(extract_guidance_info(guidance_dict=guidance_dict, attribute='SPEED_LIMIT', attribute_detail='speedLimit'))
        advisory_speed_limits.append(extract_guidance_info(guidance_dict=guidance_dict, attribute='ADVISORY_SPEED_LIMIT', attribute_detail='advisorySpeedLimit'))
        warning_signs_meanings.append(extract_guidance_info(guidance_dict=guidance_dict, attribute='WARNING_SIGN', attribute_detail='warningSign'))
        #here we get the entire dictionary containing information about these features
        overall_traffic_lights.append(extract_guidance_overall_info(guidance_dict, 'TRAFFIC_LIGHTS'))
        pedestrian_crossings.append(extract_guidance_info(guidance_dict, 'PEDESTRIAN_CROSSING', 'coordinates'))
        warning_signs.append(extract_guidance_info(guidance_dict, 'WARNING_SIGN', 'coordinates'))
    return nearest_indexes, nearest_elements, urban, motorway, link_types, averageSpeeds, speed_limits, advisory_speed_limits, roundabouts, ramps, overall_traffic_lights, pedestrian_crossings, warning_signs, warning_signs_meanings, slope_values_in_percentages, curvatures

def get_street_element_proximity(longitudes, latitudes, overall_traffic_lights, pedestrian_crossings, warning_signs, warning_signs_meanings):
    '''returns a binary status columns representing if a car is in proximity (within bounding box) of road element
        For: 
            - traffic_light_nearby_binary_indicator
            - stop_sign_nearby_binary_indicator
            - pedestrian_crossing_nearby_binary_indicator 

    keyword arguments:
    latitudes -- list of the latitude time sequence points found in the recording
    longitudes -- list of the longitude time sequence points found in the recording
    overall_traffic_lights -- list of dictionaries with the guidance attribute information for traffic lights
    pedestrian_crossings -- list of dictionaries with the guidance attribute information for pedestrian crossings
    warning_signs -- list of dictionaries with the guidance attribute information for warning signs
    warning_signs_meanings -- list representing time sequence with different warning sign meanings e.g. stop sign
    '''
    traffic_light_nearby_binary_indicator = []
    stop_sign_nearby_binary_indicator = []
    pedestrian_crossing_nearby_binary_indicator = []
    road_element_nearby_bounds = extract_bound_coordinates(overall_traffic_lights)
    pedestrian_crossing_bounds_meters = extract_bound_coordinates(pedestrian_crossings)
    stop_sign_bounds_meters = extract_stop_sign_bound_coordinates(warning_signs, warning_signs_meanings)
    for latitude, longitude in zip(latitudes, longitudes):
        meter_longitude, meter_latitude = convert_position_to_meters(longitude, latitude)
        traffic_light_nearby_binary_indicator.append(within_bounds(meter_longitude, meter_latitude, road_element_nearby_bounds))
        stop_sign_nearby_binary_indicator.append(within_bounds(meter_longitude, meter_latitude, stop_sign_bounds_meters))
        pedestrian_crossing_nearby_binary_indicator.append(within_bounds(meter_longitude, meter_latitude, pedestrian_crossing_bounds_meters))
    
    return traffic_light_nearby_binary_indicator, stop_sign_nearby_binary_indicator, pedestrian_crossing_nearby_binary_indicator

def extract_map_data(geographic_file_path: str, map_info_file_path: str, input_trip_data_folder_path: str, output_map_data_folder_path: str, signal_profile: dict):
    ''' extracts and fuses map information into the time sequence
    
    keyword arguments:
    geographic_file_path -- the path to the GeoJson file of the tile id
    map_info_file_path -- the path to the map information file
    input_trip_data_folder_path -- the path to the directory containing the trip file containing map coordinates
    output_map_data_folder_path -- the path to the directory to store the trip with the extracted map data
    '''
    #timing how long it takes
    start_time = time.time()
    homedir = str(pathlib.Path(__file__).resolve().parent.parent.parent.parent.parent)
    '''
    ATTENTION: a geojson file exported with the database inspector based on the IDs extracted with the tile_id_extraction.py script has to be saved in the _temp directory
    '''
    try:
        gdf = gpd.read_file(geographic_file_path)
    except:
        logging.error("The GeoJson file %s cannot not be read", geographic_file_path)
        sys.exit()
    end_time = time.time()
    logging.info("Reading the GeoJson took %f seconds", end_time-start_time)
    output_directory = output_map_data_folder_path
    base_directory = input_trip_data_folder_path
    for file_name in os.listdir(base_directory):
        logging.info("extracting for", file_name)
        if "cruise" not in file_name:
            df = pd.read_csv(base_directory+'/'+file_name)
            #get file location data
            latitudes = df[signal_profile["Latitude"]].to_list()
            longitudes = df[signal_profile["Longitude"]].to_list()
            #extract map information
            #timing how long it takes
            start_time = time.time()
            #
            nearest_indexes, nearest_elements, urban, motorway, link_types, averageSpeeds, speed_limits, advisory_speed_limits, roundabouts, ramps, overall_traffic_lights, pedestrian_crossings, warning_signs, warning_signs_meanings, slope_values_in_percentages, curvatures = get_map_element_instances(gdf, latitudes, longitudes)
            end_time = time.time()
            logging.info("The extraction of street attributes took %f seconds", end_time-start_time)
            #
            #timing how long it takes
            start_time = time.time()
            #
            traffic_light_nearby_binary_indicator, stop_sign_nearby_binary_indicator, pedestrian_crossing_nearby_binary_indicator = get_street_element_proximity(longitudes, latitudes, overall_traffic_lights, pedestrian_crossings, warning_signs, warning_signs_meanings)
            end_time = time.time()
            logging.info("The coordinate based extraction of street elements took %f seconds", end_time-start_time)
            #
            #combine all of the extracted map information to one numpy array and then data frame and save it to a csv
            #timing how long it takes
            start_time = time.time()
            #
            all_np = np.array([longitudes, latitudes, nearest_indexes, nearest_elements, urban, motorway, link_types, averageSpeeds, speed_limits, advisory_speed_limits, 
                stop_sign_nearby_binary_indicator, pedestrian_crossing_nearby_binary_indicator, traffic_light_nearby_binary_indicator, ramps, roundabouts, slope_values_in_percentages, curvatures]).T
            out_df = pd.DataFrame(all_np)
            out_df.to_csv(map_info_file_path, header=[signal_profile["Longitude"],signal_profile["Latitude"],'nearest_indexes', 'nearest_elements', 'urban', 'motorway', 'link_types', 'averageSpeeds',
                'speed_limits', 'advisory_speed_limits', 'stop_sign_nearby_binary_indicator', 'pedestrian_crossing_nearby_binary_indicator', 'traffic_light_nearby_binary_indicator',
                'ramps','roundabouts', 'slope_values_in_percentages', 'curvatures'])
            #fuse the data extracted from the map with the raw data
            map_df = pd.read_csv(map_info_path)
            merged = df.join(map_df, lsuffix='_raw_data', rsuffix='_map_data')
            merged = merged.drop('Unnamed: 0_raw_data', axis=1)
            merged = merged.drop('Unnamed: 0_map_data', axis=1)
            merged = merged.drop(signal_profile["Longitude"]+'_map_data', axis=1)
            merged = merged.drop(signal_profile["Latitude"]+'_map_data', axis=1)
            output_file = os.path.join(output_directory, file_name)
            merged.to_csv(output_file)
            end_time = time.time()
            logging.info("The data frame merging took %f seconds", end_time-start_time)
