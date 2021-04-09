# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 12:32:22 2021

@author: jerome
"""

# %%
#from wpg.wavefront import Wavefront
from wpg.generators import build_gauss_wavefront
from wpg.srwlib import SRWLStokes, SRWLWfr
import numpy as np
import matplotlib.pyplot as plt
import time
from wpg.wavefront import Wavefront

from math import log10, floor


# plt.style.use(['science','ieee'])

def round_sig(x, sig=2):
    if x != 0:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    else:
            return x


def getPolarisationCharacteristics(S=None,Sparam=None):
    ''''
    Return the degree of polarization, 
            degree of polarization averaged over wavefront, 
            the eccentricity, 
            orientation,  and 
            chirality 
    of the polarization ellipse.
    '''
    
    
    if S is not None:
        s0,s1,s2,s3 = getStokesParamFromStokes(S)
    else:
        if Sparam is not None:
            s0, s1, s2, s3 = Sparam
        else:
            raise ValueError("S and Sparam both None")

    # degree of polarization.  Alternatively, we could use D = Ip/(Ip + Itot) where Ip is the polarised intensity and It is the total intensity.
    D = (np.sqrt((s1**2 + s2**2 + s3**2)))/s0
    
    Davg = np.mean(D)

    # eccentricity.  e = 1 => linear polarization (no chirality) 
    e = np.sqrt(s2*np.sqrt(s1**2 + s2**2)/(1+np.sqrt(s1**2 + s2**2))) # added sqrt - jerome
    
    #making sure eccentricity agrees with s3
    # e_3 = (2*(-1 + s3**2 + np.sqrt(-1*(1 - s3**2))))/(s3**2)
    # print('Eccentricity from s3 = {}'.format(e_3))
    
    # inclination of polarization ellipse (radians)
    i = 0.5*np.arctan(s2/s1)

    # chirality, c, defined for convenience 
    s3m = np.mean(s3)
    if s3m < 0:
        c = 'ccw'
    if s3m >0:
        c = 'cw'
    if s3m ==0:
        c = None


    
    return D, Davg, e, i, c


def deg_pol(S):
    ''''
    Return the degree of polarization,  
            the eccentricity, 
            orientation, and 
            chirality 
    of the polarization ellipse.
    '''
    
    
    s0,s1,s2,s3 = S

    if type is 'linear':
        D = (np.sqrt((s1**2 + s2**2 + s3**2)))/s0
    
    return D


def getStokesParamFromStokes(stk):
    
    nTot = stk.mesh.nx*stk.mesh.ny*stk.mesh.ne
    # print('nTot before reshape: {}'.format(nTot))
    # print('arS shape : {}'.format(np.shape(stk.arS)))
    
    # plt.plot(stk.arS)
    # plt.title('arS')
    # plt.show()
    
    
    if stk.mutual > 0 :
        
        # print("TESTING IF ARRAY IS REPEATED (1-3)...")
        # s = np.reshape(stk.arS,(2,int(np.size(stk.arS)/2)))
        # print(np.nonzero(s[0]-s[1]))
        # plt.plot(s[0],label="1st")
        # plt.plot(s[1],label="2nd")
        # plt.legend()
        # plt.show()
        
        # print("TESTING IF ARRAY IS REPEATED AGAIN (neighbours)...")
        # s = np.reshape(stk.arS,(int(np.size(stk.arS)/2),2))
        # print(np.nonzero(s[:,0]-s[:,1]))
        # plt.plot(s[:,0],label="1st")
        # plt.plot(s[:,1],label="2nd")
        # plt.legend()
        # plt.show()
        
        
        # # s = np.reshape(stk.arS,(4,nTot,stk.mesh.nx,stk.mesh.ny,1))
        # s = np.reshape(stk.arS,(4,int(np.size(stk.arS)/4)))
        # print('stk.arS shape after reshape: {}'.format(np.shape(s)))
        # print("TESTING IF ARRAY IS REPEATED AGAIN (1-2)...")
        # print(np.nonzero(s[0]-s[1]))
        # plt.plot(s[0],label="1st")
        # plt.plot(s[1],label="2nd")
        # plt.legend()
        # plt.show()
        # plt.plot(s[0,:])
        # plt.title('S0')
        # plt.show()
        
        s = np.reshape(stk.arS,(4,int(np.size(stk.arS)/4)))
        # s0 = np.reshape(s[0,:,:],(stk.mesh.nx,stk.mesh.ny))
        
        S0 = np.reshape(s[0,:],(nTot,stk.mesh.nx,stk.mesh.ny,2))
        S1 = np.reshape(s[1,:],(nTot,stk.mesh.nx,stk.mesh.ny,2))
        S2 = np.reshape(s[2,:],(nTot,stk.mesh.nx,stk.mesh.ny,2))
        S3 = np.reshape(s[3,:],(nTot,stk.mesh.nx,stk.mesh.ny,2))

        
        _s0 = S0.mean(3)
        _s1 = S1.mean(3)
        _s2 = S2.mean(3)
        _s3 = S3.mean(3)
        
        # plt.imshow(_s0[int(stk.mesh.nx/2),:,:])
        # plt.title("TEST")
        # plt.show()
                
        # print("Shape of 1st average of S0: {}".format(np.shape(_s0)))
        
        s0 = abs(_s0.mean(0))
        s1 = abs(_s1.mean(0))
        s2 = abs(_s2.mean(0))
        s3 = abs(_s3.mean(0))
        
        # print("shape of 2nd averag of S0: {}".format(np.shape(s0)))
        
        # s0 = np.reshape(s[0,:,:],(stk.mesh.nx**2,stk.mesh.ny**2))
    
    else:
        s = np.reshape(stk.arS,(4,stk.mesh.nx,stk.mesh.ny))
    
        s0 = np.reshape(s[0,:,:],(stk.mesh.nx,stk.mesh.ny))
        s1 = np.reshape(s[1,:,:],(stk.mesh.nx,stk.mesh.ny))
        s2 = np.reshape(s[2,:,:],(stk.mesh.nx,stk.mesh.ny))
        s3 = np.reshape(s[3,:,:],(stk.mesh.nx,stk.mesh.ny))
    
    
    # print('S0 = {}'.format(np.shape(s0)))    
    # print('S1 = {}'.format(np.shape(s1)))    
    # print('S2 = {}'.format(np.shape(s2)))    
    # print('S3 = {}'.format(np.shape(s3)))
    
    
    
    return s0, s1, s2, s3 #, s4, s5, s6, s7

def getStokes(w, mutual=0, Fx = 1, Fy = 1):
        
    #if(isinstance(w, SRWLWfr) == False):
    #    w = SRWLWfr(w)
    stk = SRWLStokes()
    
    if mutual==0:
        stk.allocate(w.mesh.ne, w.mesh.nx, w.mesh.ny, _mutual = mutual) #numbers of points vs photon energy, horizontal and vertical positions
        stk.mesh.zStart = w.mesh.zStart #30. #longitudinal position [m] at which UR has to be calculated
        stk.mesh.eStart = w.mesh.eStart #initial photon energy [eV]
        stk.mesh.eFin = w.mesh.eFin #20000. #final photon energy [eV]
        stk.mesh.xStart = w.mesh.xStart #initial horizontal position of the collection aperture [m]
        stk.mesh.xFin = w.mesh.xFin #final horizontal position of the collection aperture [m]
        stk.mesh.yStart = w.mesh.yStart #initial vertical position of the collection aperture [m]
        stk.mesh.yFin = w.mesh.yFin #final vertical position of the collection aperture [m]  
        
    elif mutual > 0:
        
        # # Mid-point of wavefield mesh (assuming centered on 0)
        # midX = 0#(w.mesh.xFin + w.mesh.xStart)/2
        # midY = 0#(w.mesh.yFin + w.mesh.yStart)/2
        
        Sx = int(Fx*w.mesh.nx)
        Sy = int(Fy*w.mesh.ny)
        print("Sampled area (pixels): {}".format([Sx,Sy]))
        
        if Sx > 76 or Sy > 76:
            print("Error: Sampled area of wavefront is too large. Change Fx/Fy to a smaller value")
            import sys
            sys.exit()
            
        stk.allocate(w.mesh.ne, int(Fx*(w.mesh.nx)), int(Fy*(w.mesh.ny)), _mutual = mutual) #numbers of points vs photon energy, horizontal and vertical positions
        stk.mesh.zStart = w.mesh.zStart #30. #longitudinal position [m] at which UR has to be calculated
        stk.mesh.eStart = w.mesh.eStart #initial photon energy [eV]
        stk.mesh.eFin = w.mesh.eFin #20000. #final photon energy [eV]
        stk.mesh.xStart = Fx*w.mesh.xStart   #initial horizontal position of the collection aperture [m]
        stk.mesh.xFin = Fx*w.mesh.xFin #final horizontal position of the collection aperture [m]
        stk.mesh.yStart = Fy*w.mesh.yStart #initial vertical position of the collection aperture [m]
        stk.mesh.yFin = Fy*w.mesh.yFin #final vertical position of the collection aperture [m]  
        
        # print("Stokes Dimensions (xStart,xFin,yStart,yFin):")
        # print(stk.mesh.xStart)
        # print(stk.mesh.xFin)
        # print(stk.mesh.yStart)
        # print(stk.mesh.yFin)
        
        # print("Wavefield Dimensions (xStart,xFin,yStart,yFin):")
        # print(w.mesh.xStart)
        # print(w.mesh.xFin)
        # print(w.mesh.yStart)
        # print(w.mesh.yFin)
        
    """ Getting sampled area [m]"""
    Dx = Fx*w.mesh.xFin - Fx*w.mesh.xStart
    Dy = Fx*w.mesh.yFin - Fx*w.mesh.yStart
    print("Sampled range (x,y) [m]:{}".format((Dx,Dy)))
    
    w.calc_stokes(stk)
    
    # print("mutual:")
    # print(stk.mutual)
    
    return stk, Dx ,Dy

def normaliseStoke(S):
    # print("STOKE SHAPE: {}".format(np.shape(S)))
    
    s0,s1,s2,s3 = S[0], S[1], S[2], S[3]
    
    _s0 = s0/s0
    _s1 = s1/s0
    _s2 = s2/s0
    _s3 = s3/s0

    _s = np.array([[_s0.mean(),_s1.mean(),_s2.mean(),_s3.mean()]]).T
    print("Normalised Stokes vector:")
    print(_s)
    
    return _s0, _s1, _s2, _s3
    
def coherenceFromSTKS(S, Dx, Dy, pathCS = None, pathCSL = None):
    
    nTot = S.mesh.nx*S.mesh.ny*S.mesh.ne
    
    C = S.to_deg_coh()
    print("shape of coherence array: {}".format(np.shape(C)))
    
    # plt.plot(C)
    # plt.title('C')
    # plt.show()
    
    d = np.reshape(C,(nTot,S.mesh.nx,S.mesh.ny))
    print("shape of new coherence array: {}".format(np.shape(d)))
    
    dC = abs(d.mean(0))
    
    print("Shape of even newer Coherence array: {}:".format(np.shape(dC)))
    
    # print(np.nonzero(C))
    # plt.plot(C)
    # plt.show()
    
    
    Nx = int(np.squeeze(np.shape(dC[:][0])))
    Ny = int(np.squeeze(np.shape(dC[0][:])))
    
    """ Creating array of custom tick markers for plotting """
    tickAx = [round_sig(-Dx*1e6/2),round_sig(-Dx*1e6/4),0,round_sig(Dx*1e6/4),round_sig(Dx*1e6/2)]
    tickAy = [round_sig(Dy*1e6/2),round_sig(Dy*1e6/4),0,round_sig(-Dy*1e6/4),round_sig(-Dy*1e6/2)]
    
    
    plt.imshow(dC)
    plt.title("Degree of Coherence (from Stokes)")
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    plt.colorbar()
    if pathCS != None:
        print("Saving figure to path: {}".format(pathCS))
        plt.savefig(pathCS)
    plt.show()
    
    
    x0 = 0 # These are in _pixel_ coordinates!!
    y0 = 0
    x1 = int(np.max(np.shape(dC[:,0]))) # These are in _pixel_ coordinates!!
    y1 = int(np.max(np.shape(dC[0,:])))
    
    
    numx = x1 - x0 # number of points for line profile
    numy = y1 - y0
    midX = int(numx/2)
    midY = int(numy/2)
    
    X = dC[x0:x1, midY]
    Y = dC[midX, y0:y1]
    
    plt.plot(X, label="Horizontal Profile")
    plt.plot(Y, label="Vertical Profile")
    plt.legend()
    if pathCSL != None:
        print("Saving figure to path: {}".format(pathCSL))
        plt.savefig(pathCSL)
    plt.show()
    

    # cfr = w.toComplex()
    # A = cfr[x0:x1,y0:y1]
    # I = np.squeeze((abs(A.conjugate()*A)))
    
    # print("Shape of Intensity array: {}".format(np.shape(I)))
    
    # U = dC/I
    
    
    # print("Shape of U array: {}".format(np.shape(U)))
    
    # plt.imshow(U)
    # plt.title("Degree of Coherence (maybe)")
    # plt.colorbar()
    # plt.show()
    

def plotStokes(s,S,fig1='S0',fig2='S1',fig3='S2',fig4='S3', Dx=50e-6, Dy=50e-6,
               pathS0 = None, pathS1 = None, pathS2 = None, pathS3 = None,
               pathD = None, pathE = None, pathIn = None):
    
    print("Shape of S: {}".format(np.shape(S)))
    print("Shape of s: {}".format(np.shape(s)))
    
    Nx = int(np.squeeze(np.shape(s[0][:][0])))
    Ny = int(np.squeeze(np.shape(s[0][0][:])))
    
    print("Nx={}".format(Nx))
    print("Nx={}".format(Ny))
    
    """ Creating array of custom tick markers for plotting """
    tickAx = [round_sig(-Dx*1e6/2),round_sig(-Dx*1e6/4),0,round_sig(Dx*1e6/4),round_sig(Dx*1e6/2)]
    tickAy = [round_sig(Dy*1e6/2),round_sig(Dy*1e6/4),0,round_sig(-Dy*1e6/4),round_sig(-Dy*1e6/2)]
    
    
    print("plotting Stokes parameters (S0, S1, S2, S3)...")
    plt.imshow(s[0],vmin=np.min(s),vmax=np.max(s))
    plt.title(fig1)
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    if pathS0 != None:
        print("Saving S0 figure to path: {}".format(pathS0))
        plt.savefig(pathS0)
    plt.colorbar()
    plt.show()

    plt.imshow(s[1],vmin=np.min(s),vmax=np.max(s))
    plt.title(fig2)
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    if pathS1 != None:
        print("Saving S1 figure to path: {}".format(pathS1))
        plt.savefig(pathS1)
    plt.colorbar()
    plt.show()

    plt.imshow(s[2],vmin=np.min(s),vmax=np.max(s))
    plt.title(fig3)
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")   
    if pathS2 != None:
        print("Saving S2 figure to path: {}".format(pathS2))
        plt.savefig(pathS2)
    plt.colorbar()
    plt.show()

    plt.imshow(s[3],vmin=np.min(s),vmax=np.max(s))
    plt.title(fig4)
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    if pathS3 != None:
        print("Saving S3 figure to path: {}".format(pathS3))
        plt.savefig(pathS3)
    plt.colorbar()
    plt.show()
    
    
    D, Davg, e, i, c = getPolarisationCharacteristics(S=None,Sparam=s)
    
    print('Average degree of polarisation = {}'.format(Davg))
    print('Average ellipticity = {}'.format(np.mean(e)))
    print('Average inclination = {}'.format(np.mean(i)))
    print ('Chirality: {}'.format(c))
    

    plt.imshow(D)
    plt.title('Degree of polarization')
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    if pathD != None:
        print("Saving Deg of Pol figure to path: {}".format(pathD))
        plt.savefig(pathD)
    plt.colorbar()
    plt.show()
    
    plt.imshow(e)
    plt.title('Ellipticity')
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    if pathE != None:
        print("Saving ellipticity figure to path: {}".format(pathE))
        plt.savefig(pathE)
    plt.colorbar()
    plt.show()
    
    plt.imshow(i)
    plt.title('Inclination')
    plt.xticks(np.arange(0,Nx+1,Nx/4),tickAx)
    plt.yticks(np.arange(0,Ny+1,Ny/4),tickAy)
    plt.xlabel("Horizontal Position [\u03bcm]")#"(\u03bcm)")
    plt.ylabel("Vertical Position [\u03bcm]")#"(\u03bcm)")
    if pathIn != None:
        print("Saving Inclination figure to path: {}".format(pathIn))
        plt.savefig(pathIn)
    plt.colorbar()
    plt.show()

# def plotElipse(S):
    
#     from py_pol.stokes import Stokes
    
#     s0,s1,s2,s3 = S
    
    

def test():
    eMin = 10e6
    Nx = 100
    Ny = 100
    Nz = 1
    xMin = -10e-6
    xMax = 10e-6
    yMin = -10e-6
    yMax = 10e-6
    zMin = 1
    mutual = 1
    Fx = 1/2
    Fy = 1/2
    print('-----Running Test-----')
    print('-----building wavefront-----')
    w = build_gauss_wavefront(Nx,Ny,Nz,eMin/1000,xMin,xMax,yMin,yMax,1,1e-6,1e-6,1) 
    #build_gauss_wavefront()
    
    print(w)
    
    wf = Wavefront(srwl_wavefront=w)
        
    intensity = wf.get_intensity()           
    plt.imshow(intensity)
    plt.title("Intensity")
    plt.show()


    print ('-----Getting Stokes parameters-----')
    S, Dx, Dy = getStokes(w, mutual=mutual, Fx = Fx, Fy = Fy)
    s = getStokesParamFromStokes(S)
    _s = normaliseStoke(s)
    
    print("-----Plotting Stokes parameters-----")
    plotStokes(s,S, Dx=Dx, Dy=Dy)
    
    # print("-----Plotting normalised Stokes parameters-----")
    # plotStokes(_s,S,"_s0","_s1","_s2","_s3")
    
    
    print("-----Getting degree of coherence from Stokes parameters------")
    start1 = time.time()
    coherenceFromSTKS(S,Dx,Dy)
    end1 = time.time()
    print("Time taken to get degree of coherence from Stokes (s): {}".format(end1 - start1))
    
    print ('------Done------')
    
if __name__ == '__main__':
   test()


# %%
