#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 10:14:36 2021

@author: jerome
"""
# %%
def sliceThickness(lam, x, e1 = 0.1, e2 = 0.1):
    ''''
    Return the necessary thickness (z) of each slice for multi-slice
    propagation through an object (method from Li, et al)
    params:
        lam: Incident wavelength [m] (assuming monochromatic beam)
        x: Transverse pixel size [m]
        e1: Transverse sampling factor
        e2: Longitudinal sampling factor 
    Returns:
        z: minimum slice thickness in meters
    '''
    
    z = (e2*(x**2))/((e1**2)*lam)
    
    return z

# %%
def sliceNums(t, lam, x, e1=0.1, e2=0.1):
    ''''
    Return the number of slices (N) necessary for multi-slice
    propagation through an object (method from Li, et al), 
    params:
        t: Thickness of object in direction of propagation [m]
        lam: Incident wavelength [m] (assuming monochromatic beam)
        x: Transverse pixel size [m]
        e1: Transverse sampling factor
        e2: Longitudinal sampling factor 
    Returns:
        N: Minimum number of slices necessary for multi-slice propagation
    '''
    
    N = (t*lam*(e1**2))/((x**2)*e2)
    
    return N

# %%
def test():
    
    w = 6.7e-9          # incident wavelength
    dX = 2.5e-8         # transverse pixel size
    T = 90e-9           # object thickness
    
    t = sliceThickness(w,dX)
    
    print("Necessary slice thickness [m]: {}".format(t))
    
    N = sliceNums(T,w,dX)
    
    print("Number of slices needed for multi-slice: {}".format(N))
    
    print(" ")
    print("Done")

# %%
if __name__ == '__main__':
    test()
    


