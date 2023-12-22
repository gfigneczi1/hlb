#!/usr/bin/env python
"""
Sample script that uses the read_rosbag module created using
MATLAB Compiler SDK.

Refer to the MATLAB Compiler SDK documentation for more information.
"""

from __future__ import print_function
import read_rosbag
import matlab

my_read_rosbag = read_rosbag.initialize()

filePathIn = matlab.double([], size=(0, 0))
savePathIn = matlab.double([], size=(0, 0))
topicsOut = my_read_rosbag.read_rosbag(filePathIn, savePathIn)
print(topicsOut, sep='\n')

my_read_rosbag.terminate()
