#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import xarray as xr
import netCDF4 as nc
import pandas as pd
from numpy import linalg as LA
from scipy.signal import detrend

# ### Empirical Orthogonal Function (EOF) analysis
# * input: data is a two dimensional matrix (time, space) where each column is the time series of a data point
# * output: vp, var_porc, eof_p_n, exp_coef_n
# >vp: eigenvalues, vector <br>
# >var_porc: the percentage of the variance explained by each EOF mode, vector<br>
# >eof_p_n: each column is the EOF pattern (map), which are given 
# in the same order with the eigenvalues in "vp", matrix <br>
# >exp_coef_n: each column is the expansion coefficients (time series)
# associated to each EOF pattern (map), matrix<br>
# 
# * Each EOF is normalized by dividing its standard deviation. <br>
# Then, the associated time series, (the expansion coefficient), is multiplied by the same normalization factor.<br><br>
# * The EOF is calculated in an optimal way in the case that there are much more spatial points than temporal samples.<br> The algorithm is based in Preisendorfer "Principal component analysis in metereology and oceanography", page 64. 
# 

# In[1]:


def eof_n_optimize(data):
    
#     Z = data
    Z = detrend(data, axis = 0, type = 'constant')
    # n: number of temporal samples
    # p: number of points
    [n,p] = Z.shape
    
    if n>=p:
        # Accordingly to Preisendorfer's book, pages 64-66:
        # If n>=p we estimate the EOF in the State Space Setting, i.e., using Z'*Z which is p x p.
        # If n<p  we estimate the EOF in the Sample Space Setting, i.e., using Z*Z' which is n x n 

        # state space setting
        R = np.matmul(np.transpose(Z),Z)  #(p,p)
        R = R/n
        L, eof_p = LA.eigh(R) 
        vp = L
        for i in np.arange(0,L.shape[0]):
            exp_coef[:,i] = np.matmul(Z, eof_p[:,i])
    else:
        # sample space setting
        R = np.matmul(Z,np.transpose(Z))  #(n,n)
        R = R/p
        L, eof_p_f = LA.eigh(R) # L: eigen values, eof_p_f:associated eigenvectors
        vp = L
        exp_coef_b = np.empty((p,n),dtype=float)
        for i in np.arange(0,L.shape[0]):
            exp_coef_b[:,i] = np.matmul(np.transpose(Z), eof_p_f[:,i])           
        # We transform the EOF and their associated expansion coefficient (time series) in the sample 
        # space setting to those in the State space setting
        eof_p = np.empty((p,n),dtype=float)
        exp_coef = np.empty((n,n),dtype=float)
        for i in np.arange(0,L.shape[0]):
            eof_p[:,i]= exp_coef_b[:,i] / np.sqrt(np.abs(vp[i]))
            exp_coef[:,i] = np.sqrt(np.abs(vp[i])) * eof_p_f[:,i]
    #  estimate the percentage of variance explained by each EOF/eigenvector
    var_porc = (np.abs(L)/np.sum(np.abs(L)))*100    
    # normalize the EOF
    # patten is dimensionless
    # time series has the dimension/unit of data, eg. meter
    [a_1,b_1] = eof_p.shape
    [a_11,b_11] = exp_coef.shape
    eof_p_n = eof_p
    exp_coef_n=exp_coef
    
    for i in np.arange(0,b_1):
        eof_p_n[:,i] = eof_p[:,i] / np.std(eof_p[:,i])
        exp_coef_n[:,i] = exp_coef[:,i] * np.std(eof_p[:,i])
#         eof_p_n[:,i] = eof_p[:,i] * np.std(exp_coef[:,i])
#         exp_coef_n[:,i] = exp_coef[:,i] / np.std(exp_coef[:,i])


    return vp, var_porc, eof_p_n, exp_coef_n


# In[ ]:


# get_ipython().run_cell_magic('sh', '', '\njupyter nbconvert eof_n_optimize.ipynb --to python')

