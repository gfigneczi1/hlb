#!/usr/bin/env python
"""
Sample script that uses the matlab_evaluation module created using
MATLAB Compiler SDK.

Refer to the MATLAB Compiler SDK documentation for more information.
"""

from __future__ import print_function
import matlab_evaluation
import matlab

my_matlab_evaluation = matlab_evaluation.initialize()

my_matlab_evaluation.matlab_evaluation_cover(nargout=0)

my_matlab_evaluation.terminate()
