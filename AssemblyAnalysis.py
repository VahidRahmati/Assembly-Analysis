# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:39:28 2019

@author: rahmati
"""
import edit_distance as ed
import ipdb # ipdb.set_trace()
import os
#import sys
import numpy as np
import scipy.io as sio

def import_assembly(Dir_DataNames,Active_thr,movie_idx):
   
    Dir_current = os.getcwd() # pwd
    
    os.chdir(Dir_DataNames)
    #sys.path.append('D:\Dropbox\Project\PAPER\JNS 2018\Video_2018 - vr\Main analysis\Template Matching')
    from movie_names import names #import the names of the movies as "names" list variable from the movie_names.py
    os.chdir(Dir_current)
    
    name_current = names[movie_idx][11:] # the name of the current movie to analyze
    
    dir_assmblies = 'D:\Data\\2nd Project\\Data\\Data_arranged\\All\\Assemblies\\thr' + str(Active_thr) + '\\' 
    assemblies_data = sio.loadmat(dir_assmblies+ name_current + '_SGC-ASSEMBLIES.mat')
    #asmb_keys = assemblies_data.keys()
    return assemblies_data


    
def import_spike_trains(Dir_DataNames,Dir_SpikeTrains,movie_idx,nRows=13,nCols=17):
    
    Dir_current = os.getcwd() # pwd
    
    os.chdir(Dir_DataNames)
    from movie_names import names #import the names of the movies as "names" list variable from the movie_names.py
    os.chdir(Dir_current)
    
    name_current = names[movie_idx] # the name of the current movie to analyze
    eTrain_file = sio.loadmat(Dir_SpikeTrains+name_current+'.mat')
    eTrain_keys = eTrain_file.keys()
    
    a = eTrain_file['Trains_rec']['raw'][0,0]
    b = eTrain_file['Trains_rec']['raw'][0,0].reshape(nRows*nCols,80000,order='F')
    nFrames_av = eTrain_file['Trains_rec']['nFramesWithoutDrifts'][0,0][0,0] # no. of avialable frames after excluding the animal movement periods
    
    return {"nFrames_av":nFrames_av}    


#%% Read assbly detection results and compute basic measures
class AssemblyInfo(object):
    
    def __init__(self, assemblies_data, sptrains_data=None):
        # assemblies_data: the output of cell-assembly detection method using SGC for one dataset
        self.assemblies = assemblies_data['assemblies'] 
        self.activities = assemblies_data['assembly_pattern_detection']
        self.spdata = sptrains_data
        
        
    def get_ncores(self):
        return self.assemblies.shape[0] # number of core assmblies
       
      
#class Assemblies(InfoAssembly):
#
#    def __init__(self, x):
#       super(Assemblies,self).__init__(x)
#       self.core_Patchidx = self.assemblies - 1 
       
    def get_core_PatchIdx(self):
        return self.assemblies - 1
       
    def get_core_size(self):  
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
        # Sig activity aptterns of all assmblies
        return self.activities['activityPatterns'][0,0].squeeze() 

    
    def get_nsig_patterns_all(self):
        # number of all sig patterns of all aseemblie
        return self.get_sig_patterns_all().size 

        
    def get_nunits(self):
        return self.get_sig_patterns_all()[0].size

    
    def get_affinity_vec(self):
        # the vectors of cell affinities of core assemblies (each core one vector)
        return self.activities['assemblyActivityPatterns'][0,0].squeeze()

    
    def get_affinity_mat(self,nRows=13,nCols=17):
        # the spatially arranged matrices of core assemblies (each core one matrix representing one FOV)
        nCores = self.get_ncores()
        aff_list = [[]]*nCores
        for c in np.arange(nCores):
            aff_list[c] = self.get_affinity_vec()[c].reshape(nRows,nCols,order='F')       
        return aff_list
  
      
    def get_patterns_idx(self):
        # the indices of all sig "patterns" of "each" assembly (i.e. incl. also those NOT used for determining the spatial pattern of each core assembly)
        CommStrut_data = self.activities['patternSimilarityAnalysis'][0,0]['communityStructure'][0,0].copy()
        asmbPattern_idx = CommStrut_data['assignment'][0,0].squeeze().copy() # the indices of all sig patterns of each assembly  (i.e. incl. also those NOT used for determining the spatial pattern of each core assembly)
        asmbPattern_idx = asmbPattern_idx - 1 # note that the indices were importd from Matlab
        return asmbPattern_idx
    
    def get_patterns_raster(self,nRows=13,nCols=17):
        # a numpy array of nCore lists of all sig patterns of each assembly, 
        # where each list belong to one assembly and has the shape of nRows=nCells, and nCols= nSigPatternsOfThatAssembly
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
        # the indices of sig frames which used to compute the affinity matrix
        # (their spatial average gives the "priliminary" saptial pattern ...)
        return self.activities['assemblyIActivityPatterns'][0,0].squeeze() - 1 # note that the indices were imported from Matlab
    
                                                                 
    def get_count_dist(self): 
        # Pr. distribution of potential number of communities (i.e. cell assemblies)                           
        CommStrut_data = self.activities['patternSimilarityAnalysis'][0,0]['communityStructure'][0,0].copy()
        count_dist = CommStrut_data['countDistribution'][0,0].copy() 
        return count_dist

    
    def get_assemble_relfreq(self):
        # the fraction of sig patterns, for each assembly 
        freq_vec = np.zeros(self.get_ncores())
        for c in np.arange(self.get_ncores()):
            freq_vec[c] = self.get_patterns_idx()[c].size/self.get_nsig_patterns_all()   
        return freq_vec*100
        
    
    def get_assemble_freq(self):
        # the temporal frequency of each assembly, during the available recording time (freq in [1/frame])     
        nFrames_av = self.spdata['nFrames_av']
        Tfreq_vec = np.zeros(self.get_ncores()) 
        for c in np.arange(self.get_ncores()):
            Tfreq_vec[c] = self.get_patterns_idx()[c].size/nFrames_av # in [1/frame]
        return Tfreq_vec
        
    
    
        
    
    
class AssemblyMethods(AssemblyInfo):
    
    def __init__(self,assemblies_data, sptrains_data=None):
        super(AssemblyMethods,self).__init__(assemblies_data, sptrains_data=sptrains_data)
        
        
    def calc_cohess_approx(self):
        #  Cohessivenss: average coupling of each cell to the the presence of specific assembly, i.e. 
        # (Approx. approach): average of the affinities of core cells of that assembly   
        cho_approx = np.zeros(self.get_ncores())
        temp = self.get_affinity_vec() # an array of ncorr number lists
        temp = np.vstack((temp[:])) # now it is matrix with ncorr row and 221 columns
        for c in np.arange(self.get_ncores()):
            affinities = temp[c,self.get_core_PatchIdx()[c][0]].copy() #extract affinities of the core cells only
            cho_approx[c] = np.mean(affinities)
        return cho_approx, affinities   


    def calc_cohess_exact(self):
        #  Cohessivenss: average coupling of each cell to the the presence of specific assembly, i.e. 
        # (Exact approach): average of p( cell | Assembly ) over all cells of that assembly   
        cho_exact = np.zeros(self.get_ncores())
        for c in np.arange(self.get_ncores()):
            nPatterns = self.get_patterns_raster()[c].shape[1]
            Pr_nlA = np.sum(self.get_patterns_raster()[c],axis=1)/nPatterns # P(n_k|A_i) of "All" recorded individual cells
            Pr_nlA = Pr_nlA[self.get_core_PatchIdx()[c][0]] # ... of only "Core" cells
            cho_exact[c] = np.mean(Pr_nlA)
        return cho_exact, Pr_nlA 
    
    
    def calc_initiativeness(self):
        # Initiativeness (drivingness): p(assembly | cell)
        nCores = self.get_ncores()
        raster_full = np.hstack((self.get_patterns_raster()[:])) # raster matrix of all sig patterns of all assemblies
        cellActTimesTot_vec = np.sum(raster_full,axis=1)
        drivingness = [[]]*nCores
        for c in np.arange(nCores):
            coreCells = self.get_core_PatchIdx()[c][0]
            cellActTimesAsmb_vec = np.sum(self.get_patterns_raster()[c],axis=1)
            drivingness[c]=cellActTimesAsmb_vec[coreCells]/cellActTimesTot_vec[coreCells]           
        return drivingness
         
    
    def calc_spatial_overlap(self,ss_thr):
        # Compute the degree of overlap beetween the core assempbly spatial patterns, using the Dice similarity measures
        # Determine whether an assembly is the sub-assembly of another assembly, based on the spatail patterns of their cores
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

        
        return DSC, SS_main, SS
    
    
    def calc_pattern_reliability(self):
        ipdb.set_trace()
        nCores= self.get_ncores()
        ed_cores_mats = [[]]*nCores
        ed_cores_means = [[]]*nCores
        for c in np.arange(nCores):
            raster = self.get_patterns_raster()[c]
            nPatterns = raster.shape[1]
            ed_mat = np.zeros((nPatterns,nPatterns)) 
            
            ed_mat[:] = np.nan
            for i in np.arange(nPatterns):
                for j in np.arange(nPatterns):
                    ed_mat[i,j]=ed.SequenceMatcher(raster[:,i],raster[:,j]).distance()
            ed_cores_mats[c] = ed_mat
            ed_cores_means[c] = np.nanmean(ed_mat[np.tril_indices_from(ed_mat,-1)])
        return ed_cores_mats,  ed_cores_means
                     
             
     
     
    
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
