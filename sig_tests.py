# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:13:52 2019

@author: rahmati
"""
from scipy import stats

def paired(pop1,pop2,alpha_sig):
    """ paired two-tailed test """
    
    # furst check the normality if the difference pop
    diff_vec = pop2 - pop1
    _,p_normal = stats.shapiro(diff_vec) # if pval>alpha_sig, then difference population is normally distributed
    
    
    # perform a paired two-tailed test
    if p_normal>alpha_sig: # when difference is normal
        _,p_paired = stats.ttest_rel(pop1,pop2)
        test_type = 'paired two-tailed ttest'
        
    elif p_normal<alpha_sig: # ... is Not normal
        _,p_paired = stats.wilcoxon(pop1.squeeze(),pop2.squeeze(),correction=True)
        test_type = 'paired two-tailed Wilcoxon signed rank test'
        
    issig = int(p_paired<alpha_sig)
    
    return issig, p_paired, test_type
    
        
        

    
    