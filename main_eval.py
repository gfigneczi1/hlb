#This the framework cover script for all KDP-HLB related evaluations.
import sys

print(sys.path)

from preprocess import preprocess
from evaluation import Evaluate

#Preprocess layer
preprocessor = preprocess()
currentconfig = preprocessor.execute()
#Evaluation layer
Evaluate.execute(currentconfig)