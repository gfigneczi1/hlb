#################################################################################################
#                                                                                               #
#   ATTENTION: This file has to be run with a Python environment of Python 3.10 or higher       #
#   This is because it is dependent on geojson and shapely                                      #
#                                                                                               #
#################################################################################################
#os.environ['USE_PYGEOS'] = '0'
import geopandas as gpd
from shapely.geometry import LineString, Point
import numpy as np
import matplotlib.pyplot as plt
import json
import os
import pathlib
import sys
import logging
import pandas as pd
from acceleration_learning.preprocessing.map.nds_grid import create_nds_grid
from acceleration_learning.trip_file import TripFile

def extract_from_csv(file_path, signal_profile) -> list:
    '''returns list of Shapely Point objects with every point representing a gps point based on longitude and latitude of all time points in the raw data file
    
    Keyword arguments:
    file_path -- the path to the file that we want to extract the gps points from
    '''
    df = pd.read_csv(file_path)
    latitudes = df[signal_profile["Latitude"]].to_list()
    longitudes = df[signal_profile["Longitude"]].to_list()
    gps_points = []
    for latitude, longitude in zip(latitudes, longitudes):
        gps_points.append(Point(longitude, latitude))
    return gps_points

def extract_tile_id(input_trip_folder_path: str, output_geographic_data_path: str, signal_profile: dict):
    ''' extracts tile IDs for the Data Base Inspector to extract the map information into a geojson file with the trip name
    
    keyword arguments:
    input_trip_folder_path -- the path to the folder containing the trip data including map coordinates
    output_geographic_data_path -- the path to the trip geographic data folder of the project 
    '''
    observations = list()
    # All gps points are put into one observations list
    for file_name in os.listdir(input_trip_folder_path):
        trip_file_path = os.path.join(input_trip_folder_path, file_name)
        observations.extend(extract_from_csv(trip_file_path, signal_profile))
    # Makes Line String objects out of the observation points
    logging.info("The map tile id extraction found %d observations", len(observations))
    observations = LineString(observations)
    # A buffer of 30 meters around the line of observed gps points is created
    buffer = observations.buffer(3*10e-4)
    # The geometry is extracted from the lines for the observations and the buffer
    observations = gpd.GeoDataFrame(geometry=[observations])
    buffer = gpd.GeoDataFrame(geometry=[buffer])
    # The tiles are identified with 
    tiles = create_nds_grid(buffer.total_bounds[3], buffer.total_bounds[1], buffer.total_bounds[2], buffer.total_bounds[0])
    traversed_tiles = buffer.overlay(tiles)
    id_list = traversed_tiles["NDS Tile ID"].to_list()
    logging.info("The following IDs have to be input in the Custom tile selection window in the GeoJSON Export of the NDS Database Inspector:")
    logging.info(id_list)
    logging.info("The following IDs have to be input in the Custom target window in the GeoJSON Export of the NDS Database Inspector:")
    logging.info("The tile file://%s", output_geographic_data_path)
