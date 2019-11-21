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
import CondComparison as cc


# %% parameters
params = dict()
movie_idx = 0  # index of the data set (currently not applicapable)
params['Active_thr'] = 3 # importing data with activation threshold of: CoActive_thr*STD(mean[Ft])
params['Dir_DataNames'] = 'D:\Dropbox\Project\PAPER\JNS 2018\Video_2018 - vr\Main analysis\Template Matching' 
params['Dir_Assemblies'] = 'D:\Data\\2nd Project\\Data\\Data_arranged\\All\\Assemblies\\'
params['Dir_SpikeTrains'] = 'D:\Data\\2nd Project\\Data\\Data_arranged\\All\\' # dir to reconstructed event times

params['ctrlg_idx'] = [0,2,4,8,10]; # the name indices of Ctrol_g data 
params['ctrld_idx'] = [8,12,14,16,18,20]; # the name indices of Ctrol_d data

#
#spdata = aa.import_spike_trains(params['Dir_DataNames'], params['Dir_SpikeTrains'], movie_idx)
#asmb_data = aa.import_assembly(params['Dir_DataNames'], params['Dir_Assemblies'], params['Active_thr'], movie_idx)
#aa.AssemblyMethods(asmb_data,spdata).calc_transitions()


#%% import data, as instances of assembly class
instances = cc.make_instances(params)


#%% compare assembly results of Ctrol vs. Treatment

#params['method_name'] = 'get_ncores'
cc.CtrlTreat(instances,params).barplot('get_ncores','mean')

#cc.CondComparison(params)




























