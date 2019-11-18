# -*- coding: utf-8 -*-
"""
The draft script to be used later to compare neural assembly results across
different recorded movies, under different experimental condition.  

@author: Rahmati
"""
#import ipdb
import os
import sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import AssemblyAnalysis as aa


# %% parameters
movie_idx = 1  # index of the data set (currently not applicapable)
Active_thr = 3 # importing data with activation threshold of: CoActive_thr*STD(mean[Ft])
Dir_DataNames = 'D:\Dropbox\Project\PAPER\JNS 2018\Video_2018 - vr\Main analysis\Template Matching' 
Dir_SpikeTrains = 'D:\Data\\2nd Project\Data\Data_arranged\All\\' # dir to reconstructed event times


spdata = aa.import_spike_trains(Dir_DataNames, Dir_SpikeTrains, movie_idx, nFrames_orig = 80000)
asmb_data = aa.import_assembly(Dir_DataNames, Active_thr, movie_idx)


