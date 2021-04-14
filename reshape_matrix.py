#!/usr/bin/env python
# coding: utf-8


import numpy as np
import xarray as xr
import netCDF4 as nc
import pandas as pd



def reshape_matrix(grid, flag):
#     data: list of data arrary (lon, lat)
#     flag: data flags True/False, two dimensional (lon,lat)
    ind3 = len(grid)
    [ind1,ind2] = grid[0].shape
    
    new_data = np.empty((ind3,ind1*ind2), float)
    new_Matrix = np.empty((ind3,np.sum(flag)), float)
    
    for i in np.arange(0,ind3):
        new_data[i,:] = np.ravel(grid[i])
        
    if len(flag)!=ind1*ind2:
        print("Wrong dimension! Check data.")
        return
    else:
        spacenum = -1
        for i in np.arange(0,ind1*ind2):
            if(flag[i]==True):
                spacenum=spacenum+1
                new_Matrix[:,spacenum]=new_data[:,i]
        return new_Matrix
            
            


