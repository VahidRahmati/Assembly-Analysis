# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 11:51:31 2019

@author: rahmati
"""
#import ipdb
import os
import sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

#from Assembly import import_assembly as imass
#from Assembly import InfoAssembly, Assemblies, Activities

import AssemblyAnalysis as aa


# %% parameters
movie_idx = 1  # index of the data set
# importing data with activation threshold of: CoActive_thr*STD(mean[Ft])
Active_thr = 3
Dir_DataNames = 'D:\Dropbox\Project\PAPER\JNS 2018\Video_2018 - vr\Main analysis\Template Matching'
Dir_SpikeTrains = 'D:\Data\\2nd Project\Data\Data_arranged\All\\'


# import data
# Dir_current = os.getcwd() # pwd

#a = imass(Dir_DataNames,Dir_current,Active_thr,movie_idx)

spdata = aa.import_spike_trains(Dir_DataNames, Dir_SpikeTrains, movie_idx, nFrames_orig = 80000)

assdata = aa.import_assembly(Dir_DataNames, Active_thr, movie_idx)

aa.AssemblyMethods(assdata,spdata).calc_irregularity()

x = aa.AssemblyMethods(assdata, spdata).calc_irregulairty()
v = aa.AssemblyInfo(assdata, spdata).assemble_freq()
rg = aa.AssemblyMethods(assdata, spdata).calc_initiativeness()
rr = aa.AssemblyMethods(assdata).calc_cohess_approx()
cc = aa.Activities(assdata).get_patterns_raster()
l = aa.Activities(assdata).calc_cohess_approx()
t = aa.InfoAssembly(assdata).get_bad()
q = aa.Activities(assdata).assemble_rel_freq()
h = aa.Activities(assdata).get_affinity_mat()
g = aa.InfoAssembly(assdata).get_ncores()
gg = aa.Activities(assdata).get_ncores()
ggg = aa.Assemblies(assdata).get_ncores()
r = aa.Activities(assdata).get_affinity_mat()
c = aa.Assemblies(assdata).get_core_size()
d = aa.Activities(assdata).get_nsig_patterns_all()


#
#iasmb = InfoAssembly(a)
#asmbs = Assemblies(a)
#aas = Activities(a)
# print(asmbs.get_ncores())

# %% Import data
