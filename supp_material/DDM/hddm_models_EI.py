# -*- coding: utf-8 -*-
"""
Created on Fri May  1 18:06:31 2020

"""

import numpy as np
import os, fnmatch
import hddm
import pandas as pd
from patsy import dmatrix 

 


def make_model(mypath, mydata, model_name, trace_id):   
    
    print "I am using the new models"
    
       
      
    model_filename  = os.path.join(mypath, model_name, 'modelfit-md%d.model'%trace_id)
    #print model_filename
    #print model_name

    
    elif model_name == 'standard_ddm_var' :
        m = hddm.HDDM(mydata,p_outlier=0.05,include=( 'sv', 'st'))        
           
    elif model_name == 'standard_ddm_var_bias_nosv' :
        m = hddm.HDDM(mydata,p_outlier=0.05,bias=True,include=( 'st','sz'))      
    
    
    elif model_name == 'patch_v_nosv' :
        m = hddm.HDDMRegressor(mydata,{"v~costs+dacc+costs:dacc"},p_outlier=0.05,bias=True, include=( 'st','sz') )
    
    elif model_name == 'patch_a_nosv' :
        m = hddm.HDDMRegressor(mydata,{"a~costs+dacc+costs:dacc"},p_outlier=0.05,bias=True, include=( 'st','sz') )
        
    elif model_name == 'patch_a' :
        m = hddm.HDDMRegressor(mydata,{"a~costs+dacc+costs:dacc"},p_outlier=0.05,bias=True, include=( 'sv', 'st','sz') )
        
        
    elif model_name == 'vg_v_acc' :
        m = hddm.HDDMRegressor(mydata,{"v~val_diff+vmpfc+C(nb)+val_diff:vmpfc"},p_outlier=0.05, include=( 'sv', 'st') )
        
    elif model_name == 'vg_a_acc' :
        m = hddm.HDDMRegressor(mydata,{"a~val_diff+vmpfc+C(nb)+val_diff:vmpfc"},p_outlier=0.05, include=( 'sv', 'st') )
        
  
        
    return m