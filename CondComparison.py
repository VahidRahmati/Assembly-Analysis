# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:40:55 2019

@author: rahmati
"""
import ipdb
import AssemblyAnalysis as aa
import matplotlib.pyplot as plt

#def looping_master(params,method_name):
class CondComparison(object):
    
    def __init__(self,params):
        
        out_name = params['which_output']
        method_name = params['method_name']
        ctrlg_idx = params['ctrlg_idx']
        ctrld_idx = params['ctrld_idx']
        nRecords = len(ctrlg_idx+ctrld_idx)
        
        out = [[]]*nRecords
        
        for i in range(nRecords):
            
            # import event trains data
            spdata = aa.import_spike_trains(params['Dir_DataNames'], params['Dir_SpikeTrains'], i)
            
            # import assemblies results
            asmb_data = aa.import_assembly(params['Dir_DataNames'], params['Dir_Assemblies'], params['Active_thr'], i)
            
            # create an instance of the AssemblyMethods object/class
            instance = aa.AssemblyMethods(asmb_data,spdata)
            
            
            # compute the determined method/function for the instance
            result = getattr(instance,method_name)()
            
            # extract the desired output from the result variable 
            out[i] = result
            
        self.out = out
        
        
        
    def plot_mean_bars(self):
        
        1
        ipdb.set_trace()
            
               
        
        
        
        
        
        
        
        