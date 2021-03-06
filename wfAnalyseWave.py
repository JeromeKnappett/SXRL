#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 10:56:06 2021

@author: jerome
"""
import pickle
from wpg.wavefront import Wavefront
import wfCoherence as coh
import time
import wfStokes
import numpy as np
import matplotlib.pyplot as plt

# %%
def analyseWave(path,S=1,C=1,Cp=1,Fx=1/3,Fy=1/3,
                pathS0 = None, pathS1 = None, pathS2 = None, pathS3 = None,
                pathD = None, pathE = None, pathIn = None,
                pathCS = None, pathCSL = None,
                pathX = None, pathY = None, pathC = None,
                pathCL = None, pathI = None,
                pathCm = None, pathCmL = None
                ):
    ''''
    Return the Stokes parameters and degree of coherence of a pickled wavefield:
        path: File path of pickled wavefield
        S: Get Stokes parameters?                              0 - No, 1 - Yes (mutual=0) , 2 - Yes (mutual=1)
        C: Get Coherence?                                      0 - No, 1 - Yes
        Cp: Get horizontal & vertical Coherence profiles?      0 - No, 1 - Yes
        Fx: Fraction of wavefield to sample (Horizontal)
        Fy: Fraction of wavefield to sample (vertical)
    '''
    

    with open(path, 'rb') as wav:
        w = pickle.load(wav)
    
        
    wc = Wavefront(srwl_wavefront=w)
    print(" ")
    
    if S == 1:
        print ('-----Getting Stokes Parameters-----')
        start1= time.time()
        S, Dx, Dy = wfStokes.getStokes(w,0, Fx,Fy)
        s = wfStokes.getStokesParamFromStokes(S)
        _s = wfStokes.normaliseStoke(s)
        end1 = time.time()
        print('Time taken to get Stokes parameters (s): {}'.format(end1 - start1))
        
        print("(Dx,Dy):{}".format((Dx,Dy)))
        
        print(" ")
        print("-----Plotting Stokes parameters-----")
        wfStokes.plotStokes(s,S,'S0','S1','S2','S3',
                            Dx,Dy,
                            pathS0, pathS1, pathS2, pathS3,
                            pathD, pathE, pathIn)
        print(" ")
        print("-----Plotting normalised Stokes parameters-----")
        wfStokes.plotStokes(_s,S,'_s0','_s1','_s2','_s3')
    
    elif S ==2:
        print ('-----Getting Stokes Parameters-----')
        start1= time.time()
        S, Dx, Dy = wfStokes.getStokes(w,1, Fx,Fy)
        s = wfStokes.getStokesParamFromStokes(S)
        _s = wfStokes.normaliseStoke(s)
        end1 = time.time()
        print('Time taken to get Stokes parameters (s): {}'.format(end1 - start1))
        
        print("(Dx,Dy):{}".format((Dx,Dy)))
        print("S (shape): {}".format(np.shape(S)))
        
        print(" ")
        print("-----Plotting Stokes parameters-----")
        wfStokes.plotStokes(s,S,'S0','S1','S2','S3',
                            Dx,Dy,
                            pathS0, pathS1, pathS2, pathS3,
                            pathD, pathE, pathIn)
        # print(" ")
        # print("-----Plotting normalised Stokes parameters-----")
        # wfStokes.plotStokes(_s,S,'_s0','_s1','_s2','_s3')

        
        print(" ")
        print("-----Getting degree of coherence from Stokes-----")
        start2 = time.time()
        wfStokes.coherenceFromSTKS(S,Dx, Dy, pathCS, pathCSL)
        end2 = time.time()
        print('Time taken to get Coherence from Stokes parameters (s): {}'.format(end2 - start2))
        
    if C == 1:
        print(" ")
        print("-----Getting Coherence by Convolution-----")
        start3 = time.time()
        Cw, Dx, Dy = coh.Coherence(wc, Fx,Fy,
                                    pathX, pathY, pathC,
                                    pathCL, pathI) 
        end3 = time.time()
        print('Time taken to get Coherence of Cw by convolution (s): {}'.format(end3 - start3))
            
        
        print(" ")
        print("-----Plotting Coherence by Convolution")
        coh.plotCoherence(Cw,Dx,Dy,
                          pathCm, pathCmL)
    
    # if Cp ==1:
    #     start4 = time.time()
    #     Cp = coh.coherenceProfiles(wc, Fx, Fy)
    #     end4 = time.time()
    #     print('Time taken to get Coherence profiles (s): {}'.format(end4 - start4))
    
    print(" ")    
    print("Done")
    

# %%
def getPhaseGradient(wf,polarization='horizontal', axis=0):
    ''' return phase gradient along direction specified by axis for wavefront

    * ignores slices
    '''
    
    phase = np.squeeze(wf.get_phase(polarization=polarization))
    phaseGradient = np.gradient(phase,axis=axis)
            
    return phaseGradient
    
    
# %%
def poynting(lam,I,pGrad):
    
    
    k = (2*np.pi)/lam
    
    P = (1/k)*I*(pGrad+k)
    
    return P
    

# %%
def testAnalysis():
    path = 'wavefield_1.pkl'
    
    analyseWave(path, 0, 1, 1,2/5, 2/5)# 3/8,3/8)

# %%
def testPoynting():
    
    """ Loading Pickled Wavefield """
    path = '/home/jerome/Documents/MASTERS/testWave.pkl'

    with open(path, 'rb') as wav:
        w = pickle.load(wav)
    
    wf = Wavefront(srwl_wavefront=w)
    
    I = Wavefront.get_intensity(wf)
    pGrad = getPhaseGradient(wf,'total',0)
    # Ph = Wavefront.get_phase(wf)
    
    plt.imshow(I)
    plt.title("Intensity")
    plt.show()
    
    plt.imshow(pGrad)
    plt.title("Phase Gradient")
    plt.show()
    
    dX, dY = Wavefront.pixelsize(wf)
    
    # print("shape of Ph: {}".format(np.shape(Ph)))
    print("(dX,dY): {}".format((dX,dY)))
    
    
    # pGradX = np.gradient(Ph[1])#,varargs = (dX,dY),axis = (0,1))
    
    
    lamb = 6.7e-9 #wavelength of incident radiation
    
    Px = poynting(lamb,np.squeeze(I),pGrad)
    
    print("Shape of Intensity: {}".format(np.shape(I)))
    print("Shape of pGrad: {}".format(np.shape(pGrad)))
    print("shape of Poynting vector: {}".format(np.shape(Px)))
    # print("P: {}".format(Px))

    plt.imshow(Px)
    plt.title("poynting vector")      
    plt.show()
    
    poyntY = Px[:,int(np.shape(Px)[0]/2)]
    pGradY = pGrad[:,int(np.shape(pGrad)[0]/2)]
    
    plt.plot(pGradY)
    plt.title('Phase Gradient (Vertical cut)')
    plt.show()
    
    plt.plot(poyntY)
    plt.title('Poynting Vector (Vertical cut)')
    plt.show()
    
    m1 = (abs(pGrad) > 0.04) & (abs(pGrad) < 0.08)
    
    #plt.imshow(gradient[mask])
    print("shape of m1: {}".format(np.shape(m1)))
    print("shape of pGrad: {}".format(np.shape(pGrad)))
    print("shape of pGrad[m1]: {}".format(np.shape(pGrad[m1])))
    ## plt.imshow(pGrad[m1])
    # print("pGrad[m1]: {}".format(pGrad[m1]))
    # print("m1: {}".format(m1))
    # print("max pGrad: {}".format(np.max(pGrad)))
    
    m1 = m1 != True # invert mask
    pGrad[m1]=0
    plt.imshow(pGrad)
    plt.title("m=1 phase gradient")
    plt.show()

    

# %%
if __name__ == '__main__':
    testAnalysis()
    # testPoynting()
    # testComplex()


    