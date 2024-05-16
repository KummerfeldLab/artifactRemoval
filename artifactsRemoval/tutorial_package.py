from Artifact_remove import Tissue_obj
from Artifact_remove import Artifact_detect
from Artifact_remove import Artifact_remove
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy 
import os
import gzip
from statistics import mean, stdev
import copy

import pickle
import sys

from artifactsRemoval import Artifact_remove
#import artifactsRemoval.Artifact_remove
from artifactsRemoval import Tissue_obj
from artifactsRemoval import Artifact_remove
dir = "/Users/wan00232/Documents/UMNTMC-spatial/packages/detect_pkg"
test = Artifact_remove(dir = dir)

with open('data/sample.pkl', 'wb') as file:
    pickle.dump(test, file)


test.remove_border()
test.remove_edge(distance=3)
test.remove_malfunction()
test.remove_edge(distance=6)


test.review_removing()
test.simple_cleanser()

test.remove_border()
test.remove_edge(distance=3)
test.remove_malfunction()
test.remove_edge(distance=6)
test.save(test.dir)


from artifactsRemoval import Artifact_remove
dir = "/Users/wan00232/Documents/UMNTMC-spatial/Apr_2_2024/data/055_D1"
test = Artifact_remove.Artifact_remove(dir)
