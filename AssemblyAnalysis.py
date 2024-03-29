# -*- coding: utf-8 -*-
"""
The draft code for analysing neural assembly data detected from calcium imaging data.
The code mainly contains the implementation of several information-theory measures (e.g. MI, KL, Entropy),
firing statistic measures (e.g. CV2), similarity analysis (e.g. Dice similarity, edit distance), as well
as some population-coupling mesures (e.g. affinity, cohessiveness).  

The codes are currently applied to the neural assembly detection results of a single recording data (calcium 
imaging movie), and  will be later applied to the whole dataset.
 
@author: Vahid Rahmati
"""
import editdistance as ed
import ipdb # ipdb.set_trace()
import os
#import sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy import stats
import sys
import collections

def import_assembly(Dir_DataNames, Dir_Assemblies, Active_thr,movie_idx):
    """ Import the data of detected neural assmblies from MAT files""" 

    # import assembly data based on mouse ID
    sys.path.insert(0, Dir_DataNames)
    from movie_names import names    
    name_current = names[movie_idx][11:] # the name of the current movie to analyze   
    dir_assmblies = Dir_Assemblies+ 'thr' + str(Active_thr) + '\\' 
    assemblies_data = sio.loadmat(dir_assmblies+ name_current + '_SGC-ASSEMBLIES.mat')
    sigFrame_times = sio.loadmat(dir_assmblies+ name_current +'_ACTIVITY-RASTER.mat')['activity_raster_peaks']
    asmb_keys = assemblies_data.keys()    
    assemblies_data['sigFrame_times'] = sigFrame_times  
    
    
    # import sample data
#    Dir_current = os.getcwd() # pwd   
#    assemblies_data = sio.loadmat(Dir_current+'\ctrl_130618_SGC-ASSEMBLIES.mat')
#    sigFrame_times = sio.loadmat(Dir_current+'\ctrl_130618_ACTIVITY-RASTER.mat')['activity_raster_peaks']
#    assemblies_data['sigFrame_times'] = sigFrame_times
      
    return assemblies_data


    
def import_spike_trains(Dir_DataNames,Dir_SpikeTrains,movie_idx, nRows=13,nCols=17):
    """ Import the corresponding, reconstructed event trains and meta-data of each movie from MAT files"""
    
    # import event-train data based on mouse ID
    sys.path.insert(0, Dir_DataNames)
    from movie_names import names 
    name_current = names[movie_idx] # the name of the current movie to analyze
    eTrain_file = sio.loadmat(Dir_SpikeTrains+name_current+'.mat')
    eTrain_keys = eTrain_file.keys()
    drifts_mat = eTrain_file['Trains_rec']['drifts_mat'][0,0]
    nFrames_av = eTrain_file['Trains_rec']['nFramesWithoutDrifts'][0,0][0,0] # no. of avialable frames after excluding the animal movement periods
    nFrames_orig = eTrain_file['Trains_rec']['raw'][0,0].shape[2]
  
    # import sample data
#    Dir_current = os.getcwd() # pwd
#    eTrain_file = sio.loadmat(Dir_current+'\Trains_rec_ctrl_130618.mat')
#    a = eTrain_file['Trains_rec']['raw'][0,0]
#    b = eTrain_file['Trains_rec']['raw'][0,0].reshape(nRows*nCols,80000,order='F')
#    drifts_mat = eTrain_file['Trains_rec']['drifts_mat'][0,0]
#    nFrames_av = eTrain_file['Trains_rec']['nFramesWithoutDrifts'][0,0][0,0] # no. of avialable frames after excluding the animal movement periods
#    nFrames_orig = eTrain_file['Trains_rec']['raw'][0,0].shape[2]
    
    return {"nFrames_av":nFrames_av, 'drifts_mat':drifts_mat, "nFrames_orig":nFrames_orig}    


#%%
class AssemblyInfo(object):
    """ Extract the assembly detection results regarding size, frequency, spatial shape, etc. 
     In addition, compute some basic mearures sucuh as the relative, occurence frequency of each assembly.
     The main analysis are implemented in the AssemblyMethods sub-class (see below)"""

    def __init__(self, assemblies_data):
        """ assemblies_data: the output of cell-assembly detection method using SGC for one dataset"""
        
        self.assemblies = assemblies_data['assemblies'] 
        self.activities = assemblies_data['assembly_pattern_detection']
        
        
    def get_ncores(self):
        """ get the number of core assmblies"""
        
        return self.assemblies.shape[0] 
            
#class Assemblies(InfoAssembly):
#
#    def __init__(self, x):
#       super(Assemblies,self).__init__(x)
#       self.core_Patchidx = self.assemblies - 1 
       
    def get_core_PatchIdx(self):
        """ get the spatial postion of core patches of each assembly, in the FOV"""
        
        return self.assemblies - 1
       
        
    def get_core_size(self):
        """ extract the spatial size of the core of each assembly"""
        
        nCores = self.get_ncores() # number of core assmblies
        core_size = np.zeros(nCores) # size of each core assmbly ([nPaches])
        core_size[:]= np.nan
        for iC in np.arange(nCores):
            core_size[iC] = self.assemblies[iC][0].size      
        return core_size

#class Activities(InfoAssembly):
#    
#    def __init__(self,assemblies_data):
#        super(Activities,self).__init__(assemblies_data)

               
    def get_sig_patterns_all(self):
        """ get the significant activity patterns of all detected assmblies"""
        
        return self.activities['activityPatterns'][0,0].squeeze() 

    
    def get_nsig_patterns_all(self):
        """ get the number of all significant patterns of all aseemblie"""
        
        return self.get_sig_patterns_all().size 

        
    def get_nCells(self):
        """ get the number of recorded cells """ 
        
        return self.get_sig_patterns_all()[0].size

    
    def get_affinity_vec(self):
        """ get the vectors of cell affinities of core assemblies (one vector, for each core)"""
        
        return self.activities['assemblyActivityPatterns'][0,0].squeeze()

    
    def get_affinity_mat(self,nRows=13,nCols=17):
        """get the spatially arranged matrices of core assemblies (one matrix representing one FOV, for each core)"""
        
        nCores = self.get_ncores()
        aff_list = [[]]*nCores
        for c in np.arange(nCores):
            aff_list[c] = self.get_affinity_vec()[c].reshape(nRows,nCols,order='F')       
        return aff_list
  
      
    def get_patterns_idx(self):
        """ get the indices of all significant "patterns" of "each" assembly (i.e. incl. also those NOT used for determining the 
        spatial pattern of each core assembly)"""
        
        CommStrut_data = self.activities['patternSimilarityAnalysis'][0,0]['communityStructure'][0,0].copy()
        asmbPattern_idx = CommStrut_data['assignment'][0,0].squeeze().copy() # the indices of all sig patterns of each assembly  (i.e. incl. also those NOT used for determining the spatial pattern of each core assembly)
        asmbPattern_idx = asmbPattern_idx - 1 # note that the indices were importd from Matlab
        return asmbPattern_idx
    
    
    def get_patterns_raster(self,nRows=13,nCols=17):
        """ create a numpy array of nCore lists of all sig patterns of each assembly, where each list
        belongs to one assembly and has the shape of nRows=nCells, and nCols= nSigPatternsOfThatAssembly"""
        
        nCores = self.get_ncores()
        patterns_mat = [[]]*nCores
        
        for c in np.arange(nCores):
            indices = self.get_patterns_idx()[c]
            nPatterns = indices.size
            raster = np.zeros((nRows*nCols,nPatterns))
            for p in np.arange(nPatterns):
                raster[:,p] = self.get_sig_patterns_all()[p]
            patterns_mat[c] = raster.copy()            
        return patterns_mat     
        
    
    def get_affinty_patterns_idx(self):
        """ get the indices of significant frames which have been used to compute the affinity matrix
        (hint: their spatial average gives the "priliminary" saptial pattern ...)"""
        
        return self.activities['assemblyIActivityPatterns'][0,0].squeeze() - 1 # note that the indices were imported from Matlab
    
                                                                 
    def get_count_dist(self): 
        """ get probability distribution of potential number of communities (i.e. cell assemblies)"""       
                   
        CommStrut_data = self.activities['patternSimilarityAnalysis'][0,0]['communityStructure'][0,0].copy()
        count_dist = CommStrut_data['countDistribution'][0,0].copy() 
        return count_dist

    
    def get_assembly_relfreq(self):
        """ compute the fraction of significant patterns, for each assembly (relative frequency of assembly occurrence)""" 
        
        freq_vec = np.zeros(self.get_ncores())
        for c in np.arange(self.get_ncores()):
            freq_vec[c] = self.get_patterns_idx()[c].size/self.get_nsig_patterns_all()   
        return freq_vec
             

#%%            
class AssemblyMethods(AssemblyInfo):
    """ Applying measures of information-theory, similarity, sequence detection, etc, to the assess the spatiotemporal strutures 
    and dynamics of different neural assemblies detected in a single data set"""
    
    def __init__(self,assemblies_data, sptrains_data=None):
        super(AssemblyMethods,self).__init__(assemblies_data)
        
        self.nFrames_av = sptrains_data['nFrames_av']
        self.nFrames_orig = sptrains_data['nFrames_orig'] 
        
        # extract the movement-periods (drifts) of mouse during recording
        drifts_mat = sptrains_data['drifts_mat'] - 1 # indices were imported from Matlab
        drifts_mat = np.vstack((np.array([0,0]),drifts_mat)) if drifts_mat[0,0]!=0 else drifts_mat
        drifts_mat[-1,1] = self.nFrames_av if drifts_mat[-1,1] == self.nFrames_av-1 else drifts_mat[-1,1] 
        drifts_mat = np.vstack((drifts_mat,np.array([self.nFrames_orig,self.nFrames_orig]))) if drifts_mat[-1,1]!=self.nFrames_av else drifts_mat
        self.drifts_mat = drifts_mat # a matrix of the onset and ending times of drift periods 
        
        # extract the the timing or occurrence of each assembly, during reociding (based on the peak significant frames)  
        sigTimes = assemblies_data['sigFrame_times'].transpose() - 1 # indices were imported from Matlab
        sigTimes_corr = sigTimes.copy() 
        sigTimes_corr = np.double(sigTimes_corr)
        mm = 0
        for i in np.arange(drifts_mat.shape[0]):
            sigTimes_corr[drifts_mat[i,0]<sigTimes_corr] += drifts_mat[i,1] - drifts_mat[i,0] + mm            
            mm = 1
            
        self.sigTimes = sigTimes_corr    
            

    def get_labeled_times(self):
        """ assign assembly labels to the timing of significant frames of each assembly"""
        
        nCores = self.get_ncores()
        sig_times = self.sigTimes # the timing of all singinifcant frames of all assemblies
        sig_idx = self.get_patterns_idx()
        nPatterns = sig_times.size
        label_time = np.zeros((2,nPatterns))
        label_time[0,:] = sig_times
#        ipdb.set_trace()
        for i in np.arange(nCores):
            label_time[1,sig_idx[i]] = i # assmeblies are labeled as 0, 1, 2, 3, ...
      
        return label_time    


    def calc_size_corr(self):
        """ calculate the between the assembly size of the current and next one in time (all assemblies) """
        
        nCores = self.get_ncores()
        label_time = self.get_labeled_times().astype(int)
        sig_times = label_time[0,].copy()
        drifts_mat  = self.drifts_mat         
        nDrifts = drifts_mat.shape[0]
        nChunks = nDrifts - 1 # note that we have added virtual drifts to the first and end of recoding (see the importing code above)        
        count_mat = np.zeros((nCores,nCores))
        ch_size = np.zeros(nChunks)
        ch_size[:] = np.nan
        sig_patterns = self.get_sig_patterns_all().copy()
        new_label = nCores
        size_vec = []
        
        # now count the specif orders of assemblies (repetition time of all observed sequences)
        for cc in np.arange(nChunks):
            ch_start = drifts_mat[cc,1]
            ch_end = drifts_mat[cc+1,0]
            ch_indices = np.where(np.logical_and(sig_times>=ch_start, sig_times<ch_end))[0]              
            ch_label_time = label_time[:,ch_indices].copy()
            ch_label = ch_label_time[1,:]
            ch_sig = ch_label_time[0,:]           
            ch_size[cc] = ch_sig.size # the number of assemblies occuring within each chunk of time            
            if ch_sig.size>1:
                for i in np.arange(ch_sig.size):
                    size_current = sig_patterns[ch_indices[i]].sum()
                    size_vec += [size_current]
                size_vec += [np.nan]
            
        size_now = np.array(size_vec[:-1])
        size_next = np.array(size_vec[1:])
         
        NotNan_indices = np.where(np.logical_and(~np.isnan(size_now),~np.isnan(size_next)))[0]
#        ipdb.set_trace() 
        
        size_now = size_now[NotNan_indices]
        size_next = size_next[NotNan_indices]
        
        # compute the correlation (Spearman Correlation)
        corr = stats.spearmanr(size_now,size_next)
        
        return  corr[0],corr[1] #{'corr':corr[0],'pval':corr[1]} 


    def calc_transitions(self, order=1):
        """ calculate the transition matrix of in-time evolution of assemblies, where each assembly
        is considered as an state """
        
        nCores = self.get_ncores()
        label_time = self.get_labeled_times().astype(int)
        label_all = label_time[1,:].copy()
        sig_times = label_time[0,].copy()
        drifts_mat  = self.drifts_mat         
        nDrifts = drifts_mat.shape[0]
        nChunks = nDrifts - 1 # note that we have added virtual drifts to the first and end of recoding (see the importing code above)        
        count_mat = np.zeros((nCores,nCores))
        ch_size = np.zeros(nChunks)
        ch_size[:]=np.nan
        label_new = []
        
        # now count the specif orders of assemblies (repetition time of all observed sequences)
        for cc in np.arange(nChunks):
            ch_start = drifts_mat[cc,1]
            ch_end = drifts_mat[cc+1,0]              
            ch_label_time = label_time[:,np.where(np.logical_and(sig_times>=ch_start, sig_times<ch_end))[0]].copy()
            ch_label = ch_label_time[1,:]
            ch_sig = ch_label_time[0,:]           
            ch_size[cc] = ch_sig.size # the number of assemblies occuring within each chunk of time       

            if ch_sig.size>0:
                
                if ch_sig.size > 1:
                    label_new += ch_label.tolist()
                    
                for i in np.arange(ch_sig.size-order):
                    count_mat[ch_label[i],ch_label[i+order]] += 1  
        
        col_zero = np.where(count_mat.sum(axis=0)==0)
        row_zero = np.where(count_mat.sum(axis=1)==0)
        
        zero_idx = np.array([])
        if col_zero[0].size !=0:
            count_mat = np.delete(count_mat, col_zero, axis=0)
            count_mat = np.delete(count_mat, col_zero, axis=1)
            zero_idx = col_zero[0]
            
        elif row_zero[0].size!=0:
            count_mat = np.delete(count_mat, row_zero, axis=0)
            count_mat = np.delete(count_mat, row_zero, axis=1)
            zero_idx = row_zero[0]
       
        
#        counter = collections.Counter(label_new)
#        counter.values()
#        counter.keys()
#            
#        c = 0
#        count_node = np.zeros((nCores,2))
#        for k,v in counter.items():
#            count_node[c,:] = [k,v]#np.vstack(count_node,[k,v])
#            c +=1
      
        return {'count_mat':count_mat,'ch_size':ch_size,'zero_idx':zero_idx} # , 'count_node':count_node    
    
    
    def calc_assembly_seq3(self,nShuffles=5000,order=1):
        """ Compute the pvalues of assembly in-time transitions, based on transition-divergence.
        Third shuffling method: The drift timings are considered as new labels and are 
        shuffled together with assembly labels."""
        
        nCores = self.get_ncores()
        label_time = self.get_labeled_times().astype(int)
        sig_all = label_time[0,:]
        trans = self.calc_transitions(order)
        count_mat_emp = trans['count_mat']
        ch_size = trans['ch_size']
        nChunks = ch_size.size
        times_all = label_time[0,:].copy()
        labels_all = label_time[1,:].copy()
        new_labels = []
        new_times = []
        new_label = nCores # note labels start from zero
        z = 0
        for hh in np.arange(nChunks):
            ipdb.set_trace()
            if ch_size[hh]>1:
               indices =  np.arange(z, z+ch_size[hh]).astype(int)
               temp_labels = labels_all[indices]
               new_labels = np.hstack((new_labels,temp_labels,new_label))
               temp_times = times_all[indices]
               new_time = temp_times[-1] + 5  
               new_times = np.hstack((new_labels,temp_times,new_time))
               
            z += ch_size[hh]
        
        new_labels = new_labels[:-1]
        
        
        count_mat_sh = np.zeros((nCores,nCores,nShuffles))
        PrAB_sh = np.zeros_like(count_mat_sh)
        sum_mat = np.zeros_like(count_mat_emp)
        PrA_sh = np.zeros((nCores,nShuffles))
        
    
    
    def calc_assembly_seq1(self,nShuffles=5000,order=1,alpha_sig = 0.05):
        """ Compute the pvalues of assembly in-time transitions, based on transition-matrix.
        Second shuffling method: only the aseembly labeles are shuffled across chunks, and 
        drift periods remain fixed."""
        
        nCores = self.get_ncores()
        label_time = self.get_labeled_times().astype(int)
        sig_all = label_time[0,:]
        trans = self.calc_transitions(order)
        count_mat_emp = trans['count_mat']
#        nCores = count_mat_emp.shape[0]
#        count_node = trans['count_node']
        ch_size = trans['ch_size']
        zero_idx = trans['zero_idx']
        nChunks = ch_size.size
        count_mat_sh = np.zeros((nCores,nCores,nShuffles))
        PrAB_sh = np.zeros((nCores-zero_idx.size,nCores-zero_idx.size,nShuffles))
        sum_mat = np.zeros_like(count_mat_emp)
        PrA_sh = np.zeros((nShuffles,nCores-zero_idx.size))
        
        
        # for computing the count matrix in shuffled data, we don't care about the shuffles
        for s in np.arange(nShuffles):
            
            label_all = label_time[1,:].copy()
            np.random.shuffle(label_all)
            m = 0
            
            for cc in np.arange(nChunks):
                size_current = ch_size[cc]
                if size_current>0:
                    ch_label = label_all[np.arange(m,m+size_current).astype(int)]
                    m += size_current
#                    ipdb.set_trace()
                    for i in np.arange(ch_label.size-order):
                        count_mat_sh[ch_label[i],ch_label[i+order],s] += 1 
           
            
            count_sh_temp = count_mat_sh[:,:,s].copy()
            count_sh_temp =np.delete(count_sh_temp,zero_idx,axis=0)
            count_sh_temp =np.delete(count_sh_temp,zero_idx,axis=1)
            
            sum_mat += count_sh_temp>=count_mat_emp
            rows_sum_sh = np.sum(count_sh_temp,axis=1).reshape(count_sh_temp.shape[0],1) 
            temp1 = count_sh_temp/rows_sum_sh
            temp1[np.isnan(temp1)]=0
            PrAB_sh[:,:,s] = temp1
            nObserved_patterns_all = count_sh_temp.sum()
            temp2 = rows_sum_sh/nObserved_patterns_all
#            temp2 = count_node[:,1]/count_node[:,1].sum()
            temp2[np.isnan(temp2)]=0
            PrA_sh[s,:] = np.squeeze(temp2)
            
#            sum_mat += count_mat_sh[:,:,s]>=count_mat_emp
#            rows_sum_sh = np.sum(count_mat_sh[:,:,s],axis=1).reshape(count_mat_sh.shape[0],1) 
#            temp1 = count_mat_sh[:,:,s]/rows_sum_sh
#            temp1[np.isnan(temp1)]=0
#            PrAB_sh[:,:,s] = temp1
#            nObserved_patterns_all = count_mat_sh[:,:,s].sum()
#            temp2 = rows_sum_sh/nObserved_patterns_all
##            temp2 = count_node[:,1]/count_node[:,1].sum()
#            temp2[np.isnan(temp2)]=0
#            PrA_sh[s,:] = np.squeeze(temp2)
            
        # now compute pvalues
        pvals = sum_mat/nShuffles
        
        # additionally, also convert the count_mat to a Pr transition mat
        rows_sum = np.sum(count_mat_emp,axis=1).reshape(nCores-zero_idx.size,1) 
        PrAB_emp = count_mat_emp/rows_sum
        PrAB_emp[np.isnan(PrAB_emp)]=0
        
        # also compute the probabilty of observating each node (i.e. assembly: A1,A2,A3); i.e. P(s_t=A1), P(s_t=A2),P(s_t=A3), ...
        nObserved_patterns_all = count_mat_emp.sum()
        PrA_emp = rows_sum/nObserved_patterns_all
#        PrA_emp = count_node[:,1]/count_node[:,1].sum()
        PrA_emp[np.isnan(PrA_emp)] = 0
        PrA_emp = PrA_emp.reshape(1,nCores-zero_idx.size)
        ipdb.set_trace()
        nSigSeqs = pvals<alpha_sig
        nSigSeqs = nSigSeqs.sum()
        ipdb.set_trace()
        return {'pvals':pvals,'count_mat':count_mat_emp,'PrAB_sh':PrAB_sh,'PrAB_emp':PrAB_emp, 'PrA_emp':PrA_emp, 'PrA_sh':PrA_sh}
    
    
    def calc_entropy(self):
        
        PrA = self.calc_assembly_seq1()['PrA_emp']
        H_emp = - np.sum(PrA*np.log2(PrA))
#        ipdb.set_trace()
        
        return H_emp
  
    
    def calc_assembly_seq2(self,nShuffles=5000, order=1):
        """ Compute the pvalues of assembly in-time transitions, based on transition-divergence.
        First shuffling method: The drift periods are exlucded, and shuffling is performed 
        wihout considering any chunk (caveat: this will lead to a bigger number of transisitons."""
        
        nCores = self.get_ncores()
        label_time = self.get_labeled_times().astype(int)
        sig_all = label_time[0,:]
        trans = self.calc_transitions(order)
        count_mat_emp = trans['count_mat']
        count_node = trans['count_node']
        count_mat_sh = np.zeros((nCores,nCores,nShuffles))
        PrAB_sh = np.zeros_like(count_mat_sh)
        PrA_sh = np.zeros((nShuffles,nCores))
        sum_mat = np.zeros_like(count_mat_emp)
        
        # for computing the count matrix in shuffled data, we don't care about the shuffles
        for s in np.arange(nShuffles):
            label_all = label_time[1,:].copy()
            np.random.shuffle(label_all)
            for i in np.arange(sig_all.size-order):
                count_mat_sh[label_all[i],label_all[i+order],s] += 1              
            sum_mat += count_mat_sh[:,:,s]>=count_mat_emp
            rows_sum_sh = np.sum(count_mat_sh[:,:,s],axis=1).reshape(nCores,1) 
            temp1 = count_mat_sh[:,:,s]/rows_sum_sh
            temp1[np.isnan(temp1)]=0
            PrAB_sh[:,:,s] = temp1
#            nObserved_patterns_all = count_mat_sh[:,:,s].sum()
#            temp2 = rows_sum_sh/nObserved_patterns_all
            temp2 = count_node[:,1]/count_node[:,1].sum()
            temp2[np.isnan(temp2)]=0
            PrA_sh[s,:] = np.squeeze(temp2)
            
        # now compute pvalues
        pvals = sum_mat/nShuffles
        
        # additionally, also convert the count_mat to a Pr transition mat
        rows_sum = np.sum(count_mat_emp,axis=1).reshape(nCores,1) 
        PrAB_emp = count_mat_emp/rows_sum
        ipdb.set_trace()
        # also compute the probabilty of observating each node (i.e. assembly: A1,A2,A3); i.e. P(s_t=A1), P(s_t=A2),P(s_t=A3), ...
#        nObserved_patterns_all = count_mat_emp.sum()
#        PrA_emp = rows_sum/nObserved_patterns_all
        PrA_emp = count_node[:,1]/count_node[:,1].sum()
        PrA_emp = PrA_emp.reshape(1,nCores)    
        
        return {'pvals':pvals,'count_mat':count_mat_emp,'PrAB_sh':PrAB_sh,'PrAB_emp':PrAB_emp, 'PrA_emp':PrA_emp, 'PrA_sh':PrA_sh}
    
    
    def calc_KL_transitions(self,nShuffles=5000,order=1,alpha_sig = 0.05):
        """ Use KL-divergence to assess whether the state stransitions of each  node (each node separately; a node = a neural assembly)  
        is temporally more structured than a purely random process (i.e. process under uniform assumption), based on KL-divergence"""
        
#        nCores = self.get_ncores()
        seq_info = self.calc_assembly_seq1(nShuffles,order)
        PrAB_emp = seq_info['PrAB_emp']
        PrAB_sh = seq_info['PrAB_sh']
        nCores = PrAB_emp.shape[0]
        Pr_uni = np.ones_like(PrAB_emp)/nCores # the transition Pr matrix under uniform dist assumption
        
        KL_emp = np.zeros((1,nCores))
        KL_emp[:] = np.nan
        sum_mat = np.zeros_like(KL_emp)
 
        # compute the empirical KL divergence between the observed transisition Pr dist and uniform dist 
        for i in np.arange(nCores):
#            print(PrAB_emp[i,:])
            KL_emp[0,i] =  np.nansum(PrAB_emp[i,:]*np.log2(PrAB_emp[i,:]/Pr_uni[i,:]))
        
        # now use shuffled Pr mat transitions to assess whether the KL_emp is significant         
        for s in np.arange(nShuffles):
            KL_sh = np.zeros_like(KL_emp)
            KL_sh[:] = np.nan 
            for i in np.arange(nCores):
                KL_sh[0,i] =  np.nansum(PrAB_sh[i,:,s]*np.log2(PrAB_sh[i,:,s]/Pr_uni[i,:]))
            sum_mat += KL_sh>=KL_emp            
        pvals_KL = sum_mat/nShuffles   
        
        # compute the percentage of nodes with significant transitions
        nSig_nodes = pvals_KL[pvals_KL<alpha_sig].size
        nSig_nodes = nSig_nodes*100/nCores
        
        # compute the grand adverage distance of all nodes from unifrom distribution
        KL_ave = np.mean(KL_emp)
        
        return KL_emp, pvals_KL, nSig_nodes, KL_ave
    
    
    def calc_MI_transitions(self,nShuffles=5000,order=1):    
        """ Use Mutual Information to asssess whether the whole observed state transitions of all neural nodes (together) is temporally  
        more structured than that of shuffled data. MI tells us, how much knowing the current state of the process witll tell us about
        the future (next state), or vice versa. """
        
#        nCores = self.get_ncores()
        seq_info = self.calc_assembly_seq1(nShuffles,order)
        PrAB_emp = seq_info['PrAB_emp']
        PrA_emp = seq_info['PrA_emp'] 
        PrAB_sh = seq_info['PrAB_sh']  
        PrA_sh = seq_info['PrA_sh']
        nCores = PrAB_emp.shape[0]
        MI_emp = 0
        sum_mat = 0
        
        
        # (Empirical) compute the MI between current and next (future) state of whole process (i.e. all transitions of all nodes) 
        for i in np.arange(nCores):
            MI_emp += PrA_emp[0,i]*np.nansum(PrAB_emp[i,:]*np.log2(PrAB_emp[i,:]/PrA_emp))
                
        # (Suffled) compute the MI between current and next (future) state of whole process (i.e. all transitions of all nodes)   
        for s in np.arange(nShuffles):
            MI_sh = 0
            for i  in np.arange(nCores):
                MI_sh += PrA_sh[s,i]*np.nansum(PrAB_sh[i,:,s]*np.log2(PrAB_sh[i,:,s]/PrA_sh[s,:]))           
            sum_mat += MI_sh>=MI_emp       
        pval_MI = sum_mat/nShuffles        
        
        H = self.calc_entropy()
        MI_emp_norm = MI_emp/H
        
        return MI_emp, pval_MI, MI_emp_norm        
        
    
    def calc_irregularity(self):
        """ calculate the in-time irregularity of all assembly significant patterns, regardless of to whatever assembly they belong to."""
        
        sig_times = self.sigTimes # the timing of all singinifcant frames of all assemblies
        drifts_mat  = self.drifts_mat 
        nDrifts = drifts_mat.shape[0]
        nChunks = nDrifts - 1 # note that we have added virtual drifts to the first and end of recoding (see the importing code above)        
        isi_all = []
        n_isi = 0
        sum_CV2 = 0
        for cc in np.arange(nChunks):
            ch_start = drifts_mat[cc,1]
            ch_end = drifts_mat[cc+1,0]
            ch_frames = sig_times[np.where(np.logical_and(sig_times>ch_start, sig_times<ch_end))]
            
            isi_vec = np.diff(ch_frames)
            isi_all += isi_vec.tolist()
         
            
            if isi_vec.size>=2:
               for i in np.arange(isi_vec.size-1): 
#                   ipdb.set_trace()
                   sum_CV2 += 2*np.absolute(isi_vec[i]-isi_vec[i+1])/(isi_vec[i]+isi_vec[i+1])
               n_isi += isi_vec.size
           

        CV2 = sum_CV2/(n_isi-1)
        isi_all = np.array(isi_all)
        CV = np.nanstd(isi_all)/np.nanmean(isi_all)        
        
        return CV2, CV, isi_all
                    
    
    def calc_assembly_freq(self,fs=74.14):
        """ calculate the temporal frequency of each assembly, during the available recording time (freq in [1/frame])"""
        
        Tfreq_vec = np.zeros(self.get_ncores()) 
        for c in np.arange(self.get_ncores()):
            Tfreq_vec[c] = self.get_patterns_idx()[c].size/self.nFrames_av # in [1/frame]
        return Tfreq_vec*60*fs # in [min]           
    
    
    def calc_cohess_approx(self):
        """  calculate the cohessivenss: average coupling of each cell to the the presence of specific assembly; 
        (Approx. approach): average of the affinities of core cells of that assembly"""   
        
        cho_approx = np.zeros(self.get_ncores())
        temp = self.get_affinity_vec() # an array of ncorr number lists
        temp = np.vstack((temp[:])) # now it is matrix with ncorr row and 221 columns
        for c in np.arange(self.get_ncores()):
            affinities = temp[c,self.get_core_PatchIdx()[c][0]].copy() #extract affinities of the core cells only
            cho_approx[c] = np.mean(affinities)
            
        cho_ave = np.mean(cho_approx)
        
        return cho_approx, affinities, cho_ave   


    def calc_cohess_exact(self):
        """ calculate the cohessivenss: average coupling of each cell to the the presence of specific assembly; 
        (Exact approach): average of p( cell | Assembly ) over all cells of that assembly"""
        
        cho_exact = np.zeros(self.get_ncores())
        for c in np.arange(self.get_ncores()):
            nPatterns = self.get_patterns_raster()[c].shape[1]
            Pr_nlA = np.sum(self.get_patterns_raster()[c],axis=1)/nPatterns # P(n_k|A_i) of "All" recorded individual cells
            Pr_nlA = Pr_nlA[self.get_core_PatchIdx()[c][0]] # ... of only "Core" cells
            cho_exact[c] = np.mean(Pr_nlA)
        
        cho_ave = np.mean(cho_exact)
        
        return cho_exact, Pr_nlA, cho_ave 
    
    
    def calc_initiativeness(self):
        """ calculate the initiativeness: p(assembly | cell)"""
        
        nCores = self.get_ncores()
        raster_full = np.hstack((self.get_patterns_raster()[:])) # raster matrix of all sig patterns of all assemblies
        cellActTimesTot_vec = np.sum(raster_full,axis=1)
        Pr_Aln = [[]]*nCores
        Pr_Aln_ave = np.zeros(self.get_ncores())
        
        for c in np.arange(nCores):
            coreCells = self.get_core_PatchIdx()[c][0]
            cellActTimesAsmb_vec = np.sum(self.get_patterns_raster()[c],axis=1)
            Pr_Aln[c]=cellActTimesAsmb_vec[coreCells]/cellActTimesTot_vec[coreCells]   
            Pr_Aln_ave[c] = np.mean(Pr_Aln[c])
        
        ini_ave = np.nanmean(Pr_Aln_ave)
            
        return Pr_Aln_ave, Pr_Aln, ini_ave
         
    
    def calc_spatial_overlap(self,ss_thr=80):
        """ calculate the degree of overlap beetween the core assempbly spatial patterns, using the Dice similarity measures.
        Also, determine whether an assembly is a sub-assembly of another assembly, based on the spatail patterns of their cores"""
        
        nCores = self.get_ncores()
        DSC = np.zeros((nCores,nCores))
        SS = np.zeros((nCores,nCores))
        for i in np.arange(nCores):
            for j in np.arange(nCores):
                
                vec1 = self.get_core_PatchIdx()[i][0]
                vec2 = self.get_core_PatchIdx()[j][0]

                # Calculate the spatial overlap using Dice similarity
                DSC[i,j]= 2*np.intersect1d(vec1,vec2).size/(vec1.size + vec2.size)       
                
                # Determine the sub-assembliness (SS)
                SS[i,j] = np.intersect1d(vec1,vec2).size*100/vec1.size                
        
        # now determine whether e.g. core of assembly A is a sub-assmbly of that of B        
        SS_main = np.zeros((nCores, nCores)) # binary
        for i in np.arange(nCores):
            for j in np.arange(nCores):
                if SS[i,j]>=ss_thr and SS[j,i]<ss_thr:
                    SS_main[i,j] = 1
                    
#        ipdb.set_trace()
        DSC_ave = np.mean(DSC[np.tril_indices(DSC.shape[0],k=-1)])
        overlap_ave = np.mean(SS[np.tril_indices(SS.shape[0],k=-1)])
        nFull_overlap = SS_main.sum()
        
        return DSC_ave, overlap_ave, nFull_overlap
    
    
    def calc_pattern_reliability(self):
        """ calculate the average Edist-distance between the spatial pattern of each core assembly and each of its significant patterns"""  
        
        nCores= self.get_ncores()
        ed_cores_mats = [[]]*nCores
        ed_cores_means = [[]]*nCores
        core_pidx = self.get_core_PatchIdx()
        
        for c in np.arange(nCores):
            raster = self.get_patterns_raster()[c]
            nPatterns = raster.shape[1]
            ed_mat = np.zeros(nPatterns) 
            ed_mat[:] = np.nan
           
            core_binary = np.zeros(raster.shape[0])
            core_binary[core_pidx[c][0][0]] = 1
            
            for i in np.arange(nPatterns):             
                ed_mat[i]=ed.distance(core_binary.tolist(),raster[:,i].tolist())
  
            ed_cores_mats[c] = ed_mat
            ed_cores_means[c] = np.nanmean(ed_mat)
            
        ed_grandmean = np.nanmean(ed_cores_means)
        
        return ed_cores_mats, ed_cores_means, ed_grandmean    
    
    
#%% draft code    
    
    
#    os.chdir(Dir_DataNames)
#    from movie_names import names #import the names of the movies as "names" list variable from the movie_names.py
#    os.chdir(Dir_current)    
                     

#    def calc_pattern_reliability(self):
#        # calculate the average Edist-distance between the core assembly pattern and each of its sig patterns  
#        nCores= self.get_ncores()
#        ed_cores_mats = [[]]*nCores
#        ed_cores_means = [[]]*nCores
#        for c in np.arange(nCores):
#            raster = self.get_patterns_raster()[c]
#            nPatterns = raster.shape[1]
#            ed_mat = np.zeros((nPatterns,nPatterns)) 
#            ed_mat[:] = np.nan
#            for ii in np.arange(1,nPatterns,1):
#                for jj in np.arange(ii):
#                    ed_mat[ii,jj]=ed.SequenceMatcher(raster[:,ii].tolist(),raster[:,jj].tolist()).distance()
#            ed_cores_mats[c] = ed_mat
#            ed_cores_means[c] = np.nanmean(ed_mat[np.tril_indices_from(ed_mat,-1)])
#        return ed_cores_mats,  ed_cores_means     
     
    
#%% class of info_assembly
#class info_assembly(object):
#    
#    def __init__(self,x):
#        self.assemblies_data = x # assemblies_data: the output of cell-assembly detection method using SGC for one dataset
#        
#    
#    #%% Extract the patch indices of cores of assmblies    
#    def assemblies(self):        
#        core_Patchidx = self.assemblies_data['assemblies'] # patch indices of each core assmbly
#        core_Patchidx = core_Patchidx  - 1 # note that the indices were imported from Matlab
#        nCores = core_Patchidx.size # number of core assmblies
#        
#        core_size = np.empty((nCores)) # size of each core assmbly ([nPaches])
#        for iC in np.arange(nCores):
#            core_size[iC] = core_Patchidx[iC][0].size
#        
#        return {'core_Patchidx':core_Patchidx, 'nCores':nCores, 'core_size':nCores} 
#    
#    
#    #%% Extract activity patterns and strucutres 
#    def activities(self):
#        coreActivityPatterns_data = self.assemblies_data['assembly_pattern_detection']
#        
#        sig_patterns_all = coreActivityPatterns_data['activityPatterns'][0,0].squeeze() # Sig activity aptterns of all assmblies
#        nSig_patterns_all = sig_patterns_all.size # number of all sig patterns of all aseemblies
#        
#        nPatches = sig_patterns_all[0].size # number of all patches or cells recorded
#        
#        affinity_mat = coreActivityPatterns_data['assemblyActivityPatterns'][0,0].squeeze() # the affinity matrix of each core assembly
#        affintyPatterns_idx = coreActivityPatterns_data['assemblyIActivityPatterns'][0,0].squeeze() # the indices of sig frames which used to compute the affinity matrix
#                                                                                                    # (their spatial average gives the "priliminary" saptial pattern ...)
#        affintyPatterns_idx = affintyPatterns_idx - 1 # note that the indices were imported from Matlab
#        
#        CommStrut_data = coreActivityPatterns_data['patternSimilarityAnalysis'][0,0]['communityStructure'][0,0]
#        count_dist = CommStrut_data['countDistribution'][0,0] # Pr. distribution of potential number of communities (i.e. cell assemblies)
#        asmbPattern_idx = CommStrut_data['assignment'][0,0].squeeze() # the indices of all sig patterns of each assembly  (i.e. incl. also those NOT used for determining the spatial pattern of each core assembly)
#        asmbPattern_idx = asmbPattern_idx - 1 # note that the indices were importd from Matlab
#        
#        #eTrain2D_current = eTrain_file['Trains_zeroing_2D'].transpose()
#        
#        return {'sig_patterns_all':sig_patterns_all,'nSig_patterns_all':nSig_patterns_all,'asmbPattern_idx':asmbPattern_idx,
#                'nPatches':nPatches, 'affinity_mat':affinity_mat,'affinity_mat':affinity_mat, 'count_dist':count_dist}
#    
#
#    
#    
    
    
    
    
    
    
    
    
    
    
    

# =============================================================================
# def import_assembly(Dir_DataNames,Dir_current,Active_thr,movie_idx):
#    
#     os.chdir(Dir_DataNames)
#     #sys.path.append('D:\Dropbox\Project\PAPER\JNS 2018\Video_2018 - vr\Main analysis\Template Matching')
#     from movie_names import names #import the names of the movies as "names" list variable from the movie_names.py
#     os.chdir(Dir_current)
#     
#     name_current = names[movie_idx][11:] # the name of the current movie to analyze
#     
#     dir_assmblies = 'D:\Data\\2nd Project\\Data\\Data_arranged\\All\\Assemblies\\thr' + str(Active_thr) + '\\' 
#     assemblies_data = sio.loadmat(dir_assmblies+ name_current + '_SGC-ASSEMBLIES.mat')
#     asmb_keys = assemblies_data.keys()
#     
#     
#     #%% Extract Core Assmblies data
#     core_Patchidx = assemblies_data['assemblies'] # patch indices of each core assmbly
#     core_Patchidx = core_Patchidx  - 1 # note that the indices were imported from Matlab
#     nCores = core_Patchidx.size # number of core assmblies
#     
#     core_size = np.empty((nCores)) # size of each core assmbly ([nPaches])
#     for iC in np.arange(nCores):
#         core_size[iC] = core_Patchidx[iC][0].size
#     
#     
#     #%% Extract activity patterns and strucutres 
#     coreActivityPatterns_data = assemblies_data['assembly_pattern_detection']
#     
#     sig_patterns_all = coreActivityPatterns_data['activityPatterns'][0,0].squeeze() # Sig activity aptterns of all assmblies
#     nSig_patterns_all = sig_patterns_all.size # number of all sig patterns of all aseemblies
#     
#     nPatches = sig_patterns_all[0].size # number of all patches or cells recorded
#     
#     affinity_mat = coreActivityPatterns_data['assemblyActivityPatterns'][0,0].squeeze() # the affinity matrix of each core assembly
#     affintyPatterns_idx = coreActivityPatterns_data['assemblyIActivityPatterns'][0,0].squeeze() # the indices of sig frames which used to compute the affinity matrix
#                                                                                                 # (their spatial average gives the "priliminary" saptial pattern ...)
#     affintyPatterns_idx = affintyPatterns_idx - 1 # note that the indices were imported from Matlab
#     
#     CommStrut_data = coreActivityPatterns_data['patternSimilarityAnalysis'][0,0]['communityStructure'][0,0]
#     count_dist = CommStrut_data['countDistribution'][0,0] # Pr. distribution of potential number of communities (i.e. cell assemblies)
#     asmbPattern_idx = CommStrut_data['assignment'][0,0].squeeze() # the indices of all sig patterns of each assembly  (i.e. incl. also those NOT used for determining the spatial pattern of each core assembly)
#     asmbPattern_idx = asmbPattern_idx - 1 # note that the indices were importd from Matlab
#     
#     #eTrain2D_current = eTrain_file['Trains_zeroing_2D'].transpose()
#     
#     return {'nCores':nCores, 'core_size':nCores, 'sig_patterns_all':sig_patterns_all,'nSig_patterns_all':nSig_patterns_all,
#             'nPatches':nPatches, 'affinity_mat':affinity_mat,'affinity_mat':affinity_mat, 'count_dist':count_dist,
#             'asmbPattern_idx':asmbPattern_idx}
#     
# =============================================================================
