# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 16:40:55 2019

@author: rahmati
"""
import ipdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import shapiro

import AssemblyAnalysis as aa
import sig_tests as stest


def make_instances(params):
    
    ctrlg_idx = params['ctrlg_idx']
    ctrld_idx = params['ctrld_idx']
    nRecords = 2*(len(ctrlg_idx)+len(ctrld_idx))
        
    instances = [[]]*nRecords
        
    for i in range(nRecords):
            
        # import event trains data
        spdata = aa.import_spike_trains(params['Dir_DataNames'], params['Dir_SpikeTrains'], i)
            
        # import assemblies results
        asmb_data = aa.import_assembly(params['Dir_DataNames'], params['Dir_Assemblies'], params['Active_thr'], i)
            
        # create an instance of the AssemblyMethods object/class
        instances[i] = aa.AssemblyMethods(asmb_data,spdata)    
#        ipdb.set_trace()
    return instances
   

def make_panda(result,stat_name):
    
    # bring results to panda format
    stat_all = np.concatenate(
            (result['ctrld'], result['diaz'], result['ctrlg'], result['gabaz']), axis=0)
    
    # convert these data to panda structure
    pd_stat = pd.DataFrame({stat_name: stat_all[:, 0]})
    
    # assign label of condition to each row
    nCtrld = result['ctrld'].size
    nDiaz = result['diaz'].size
    nCtrlg = result['ctrlg'].size
    nGabaz = result['gabaz'].size
    Ctrld_label = np.repeat(["Control(D)"], nCtrld).reshape(-1, 1)
    Diaz_label = np.repeat(["Diazepam"], nDiaz).reshape(-1, 1)
    Ctrlg_label = np.repeat(["Control(G)"], nCtrlg).reshape(-1, 1)
    Gabaz_label = np.repeat(["Gabazine"], nGabaz).reshape(-1, 1)
    Condition_labels = np.concatenate(
        (Ctrld_label, Diaz_label, Ctrlg_label, Gabaz_label))
    
    
    # assign label of Ctrl vs. Treatment to each row
    versusCD_label = np.repeat(["Diazepam vs. Control"],
                               nCtrld+nDiaz).reshape(-1, 1)
    versusCG_label = np.repeat(["Gabazine vs. Control"],
                               nCtrlg+nGabaz).reshape(-1, 1)
    Versus_labels = np.concatenate((versusCD_label, versusCG_label))
    
    
    # add these label to our panda strucutre defined some lines above
    pd_stat['Condition'] = Condition_labels
    pd_stat['Versus'] = Versus_labels    
    
    return pd_stat



    
#def looping_master(params,method_name):
class CtrlTreat(object):
    def __init__(self,instances,params):
        self.instances = instances
        self.params= params
        
        
    def master_loop(self,method_name,out_idx=None):
        
#        method_name = params['method_name']
        ctrlg_idx = self.params['ctrlg_idx']
        ctrld_idx = self.params['ctrld_idx']
        
        out_ctrld = [[]]*len(ctrld_idx)
        out_diaz = [[]]*len(ctrld_idx)
        out_ctrlg = [[]]*len(ctrlg_idx)
        out_gabaz = [[]]*len(ctrlg_idx)
        idx = 0
        
        if out_idx is None:
                
            for i in ctrld_idx:           
                # compute the determined method/function for the instance
                out_ctrld[idx] = [getattr(self.instances[i],method_name)()]
                out_diaz[idx] = [getattr(self.instances[i+1],method_name)()]
                idx += 1
            
            idx = 0
            for i in ctrlg_idx:            
                # compute the determined method/function for the instance
                out_ctrlg[idx] = [getattr(self.instances[i],method_name)()]
                out_gabaz[idx] = [getattr(self.instances[i+1],method_name)()]  
                idx += 1
        
        elif out_idx is not None:
            
            for i in ctrld_idx:           
                # compute the determined method/function for the instance
                out_ctrld[idx] = [getattr(self.instances[i],method_name)()[out_idx]]
                out_diaz[idx] = [getattr(self.instances[i+1],method_name)()[out_idx]]
                idx += 1
            
            idx = 0
            for i in ctrlg_idx:         
#                ipdb.set_trace()
                # compute the determined method/function for the instance
                out_ctrlg[idx] = [getattr(self.instances[i],method_name)()[out_idx]]
                out_gabaz[idx] = [getattr(self.instances[i+1],method_name)()[out_idx]]  
                idx += 1
                
        
        return {'ctrld':out_ctrld, 'diaz':out_diaz, 'ctrlg':out_ctrlg, 'gabaz':out_gabaz}
        
    
    
    def cmpt_stat(self,out,stat_name='mean'):
#       out =  self.master_loop(method_name,out_idx)
       
#       ipdb.set_trace()
       stat_dict = {}
       for k,v in out.items():
           temp_stat = np.zeros((len(v),1))
           for j in range(len(v)):
               temp_stat[j] = getattr(np.array(v[j]),stat_name)() # e.g. mean over assmblies of each recording, for each condition (e.g. Ctrld)
           stat_dict[k]=temp_stat      
#           stat_dict2[k] = np.array([getattr(np.array(v[j]),stat_name)().reshape(-1,1) for j in range(len(v))])
           
       return stat_dict
       
        
#    def plot_bar(self,result, plot_type='box', alpha_sig = 0.05):
    def plot_bar(self,result, plot_type='box', alpha_sig = 0.05,title = 'measure', ylabel = 'measures', ylim = 0 ):  
        
        # check whether there is a signigicant difference between Ctrl vs. Treatment
        issig_d, pval_d, test_type_d = stest.paired(result['ctrld'],result['diaz'], alpha_sig)
        issig_g, pval_g, test_type_g = stest.paired(result['ctrlg'],result['gabaz'], alpha_sig)
        
        print('Diaz: pval=%.4f, test=%s' % (pval_d, test_type_d))
        print('Gabaz: pval=%.4f, test=%s' % (pval_g, test_type_g))
        
        
        # take the results to panda format
        pd_stat = make_panda(result,'measure')
        
        
        clrs = ['Skyblue', 'Red', 'SkyBlue', 'Red']

        # box plot
        sns.boxplot(data=pd_stat,  x="Condition", y= 'measure',
                    palette=clrs, boxprops={'facecolor': 'None', 'linewidth': 2})
        
        # Swarm plot
        ax_sw = sns.swarmplot(data=pd_stat,  x="Condition",
                              y='measure', size=10, zorder=0, dodge=True)
        
        # Add statistical annotations
        x1, x2, x3, x4 = 0, 1, 2, 3
        max_val = pd_stat['measure'].max() 
        max_val += max_val*10/100;
        
        y, h, color = max_val, max_val/30, 'k'
        
        
        # for Ctrld vs. Diaz
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, color=color)
        svs = ['*','center',30] if issig_d==1 else ['n.s.', 'bottom',15] 
        plt.text(0.5*(x1+x2), y+h, s=svs[0], ha='center',
                 va=svs[1], size=svs[2], color=color)
        
        # for Ctrlg vs. Gabaz
        plt.plot([x3, x3, x4, x4], [y, y+h, y+h, y], lw=2, color=color)
        svs = ['*','center',30] if issig_g==1 else ['n.s.', 'bottom', 15]
        plt.text(0.5*(x3+x4), y+h, s=svs[0], ha='center',
                 va=svs[1], size=svs[2], color=color)
        
#        ax_sw.axes.set_ylim(0, pd_stat[stat_name].max()+30)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.xticks(rotation=45)
        
        y_lim = plt.gca().get_ylim()
        y_max = 1.08*y_lim[1]  
        plt.ylim((y_lim[0],y_max))
        
        if ylim != 0:
            plt.ylim(ylim)
            
        ipdb.set_trace()
        
#        plt.ylim((0,12))
        plt.show()
        
        
    def plot_CV2(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Local irregularity'
        keyarg_plot['ylabel'] = 'CV2'
        
        # extract the desired measure for each mouse of each condition
        out =  self.master_loop('calc_irregularity',0)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)
        
        
    def plot_CV(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Global irregularity'
        keyarg_plot['ylabel'] = 'CV'
        
        # extract the desired measure for each mouse of each condition
        out =  self.master_loop('calc_irregularity',1)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot) 
        
    
    def plot_meanIAI(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Inter-Assembly Intervals'
        keyarg_plot['ylabel'] = 'IAI [frame]'
        
        # extract the desired measure for each mouse of each condition
        out =  self.master_loop('calc_irregularity',2)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot) 
            

    def plot_relfreq(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Relative frequency'
        keyarg_plot['ylabel'] = 'Occurrence prob.'
        
        # extract the desired measure for each mouse of each condition
        out =  self.master_loop('get_assembly_relfreq')
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)         


    def plot_freq(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'In-time frequency'
        keyarg_plot['ylabel'] = 'Occurrence frequency [1/min]'
        
        # extract the desired measure for each mouse of each condition
        out =  self.master_loop('calc_assembly_freq')
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)         


    def plot_size_corr(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Assembly size correlation'
        keyarg_plot['ylabel'] = 'Correlation'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_size_corr',1)
        
#        ipdb.set_trace()
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)         


    def plot_MI(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Temporal structure'
        keyarg_plot['ylabel'] = 'Mutula information [bit]'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_MI_transitions',0)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)    
        
        
    def plot_H(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Temporal structure'
        keyarg_plot['ylabel'] = 'Entropy [bit]'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_entropy')
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)    


    def plot_KL_SigNodes_rel(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Significant node transitions'
        keyarg_plot['ylabel'] = 'precntage of nSig nodes'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_KL_transitions',2)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)   
 

    def plot_KL_ave(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Node transitions'
        keyarg_plot['ylabel'] = 'Ave. KL distance'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_KL_transitions',3)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)
        
 
    def plot_cohess_exact(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Cohessivenss (Exact)'
        keyarg_plot['ylabel'] = 'Ave. cohessiveness'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_cohess_exact',2)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)
        

    def plot_cohess_approx(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Cohessivenss (Approximate)'
        keyarg_plot['ylabel'] = 'Ave. cohessiveness'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_cohess_approx',2)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)

 
    def plot_init(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Initiativeness'
        keyarg_plot['ylabel'] = 'Ave. initiativeness'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_initiativeness',2)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)
 

    def plot_overlap(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Assembly spatial overlap'
        keyarg_plot['ylabel'] = 'Ave. overlap'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_spatial_overlap',0)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot)       
        
        
    def plot_reliability(self, keyarg_compt=None,keyarg_plot=None,**keyargs):    
        
        keyarg_compt = {} if keyarg_compt is None else keyarg_compt
        keyarg_plot = {} if keyarg_plot is None else keyarg_plot            
        
        # add title and label to the keyarg_plot
        keyarg_plot['title'] = 'Pattern reliability'
        keyarg_plot['ylabel'] = 'Ave. of edit distances'
        
        # extract the desired measure for each mouse of each condition
        out = self.master_loop('calc_pattern_reliability',2)
        
        # compute the deisred descriptive measure for each condition
        result = self.cmpt_stat(out,**keyarg_compt)
        
        # plot the resutls
        self.plot_bar(result, **keyarg_plot) 
        
        
#    def plot_nCores(self, )
        
#     def plot_CV2(self, keyarg_compt,keyarg_plot,**keyargs):    
#        ipdb.set_trace()
#        
#        
#        # extract the desired measure for each mouse of each condition
#        out =  self.master_loop('calc_irregularity',0)
#        
#        # compute the deisred descriptive measure for each condition
#        result = self.cmpt_stat(out,**keyarg)
#        
#        # plot the resutls
#        self.plot_bar(result, keyarg)       
        
        
        
        
        
        
        