#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 19:41:40 2018

@author: tylerclark
"""
from evolve.generate import explore

    
    
parms={'num_seasons':100,
       'beta':2,
       'epsilon':None,
       'alpha':None,
       'p':0.4,
       'b':0.08,
       'n':6,
       'H':[0.33,0.34,0.33],
       'x':0.5,
       'rand_beta':[False,1],          
       'plot_dist_season':False,
       'plot_mean_var':True,
       'plot_dist_end':False,
       }

var1={'var':'p','min':0.2,'max':0.2,'steps':10}
var2={'var':'x','min':0.2,'max':1,'steps':1}


explore(parms,var1,var2)

