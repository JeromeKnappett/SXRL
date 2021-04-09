#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 10:40:11 2021

@author: jerome
"""
import numpy as np
import matplotlib.pyplot as plt
from wpg.generators import build_gauss_wavefront
from wpg.wavefront import Wavefront



from PIL import Image
from PIL import ImageDraw
import cv2

# plt.style.reload_library()
# plt.style.use(['science','no-latex'])#,'ieee'])

def gratingEfficiencyDEL(m,beta,delta,d,k,b,theta = 0):
    """ DOES NOT WORK - NEED TO CHECK """
    
    ''''
    Return the diffraction efficiency of X-rays incident on an arbitrary grating
    with rectangular wire cross sections (method from Delvaille et al. (1980)):
        m: order of diffracted beam
        beta: imaginary part of refractive index
        delta: (1-real part) of refractive index
        d: grating periodicity
        k: photon energy in units of (hbar*c)
        b: thickness of mask
        theta: incident angle
        
    '''
    a = d # square grating
    
    phi = (2*np.pi*m)/(d*np.cos(theta))
    print("phi: {}".format(phi))
    

    X_0 = 0
    X_1 = a*np.cos(theta) - b*np.sin(theta)
    X_2 = a*np.cos(theta)
    X_3 = d*np.cos(theta) - b*np.sin(theta)
    X_4 = d*np.cos(theta)
    
    Z_0 = 0
    Z_1 = 0
    Z_2 = b/np.cos(theta)
    Z_3 = b/np.cos(theta)
    
    C_0 = 0
    C_1 = 1/(np.sin(theta)*np.cos(theta))
    C_2 = 0
    C_3 = -1/(np.sin(theta)*np.cos(theta))
    
    R_0 = (1/k)*(((C_0*k*beta)**2 + (k*C_0*delta + phi)**2)**(1/2))
    R_1 = (1/k)*(((C_1*k*beta)**2 + (k*C_1*delta + phi)**2)**(1/2))
    R_2 = (1/k)*(((C_2*k*beta)**2 + (k*C_2*delta + phi)**2)**(1/2))
    R_3 = (1/k)*(((C_3*k*beta)**2 + (k*C_3*delta + phi)**2)**(1/2))

    F_0 = (1/k*R_0)*np.exp(-1j*np.cos((beta*C_0)/R_0)+1j*k*(Z_0 - C_0 * X_0)*(delta + 1j*beta))*(np.exp(1j*X_1*(k*C_0*(delta+1j*beta)+phi))-(np.exp(1j*X_0*(k*C_0*(delta+1j*beta)+phi))))
    F_1 = (1/k*R_1)*np.exp(-1j*np.cos((beta*C_1)/R_1)+1j*k*(Z_1 - C_1 * X_1)*(delta + 1j*beta))*(np.exp(1j*X_2*(k*C_1*(delta+1j*beta)+phi))-(np.exp(1j*X_1*(k*C_1*(delta+1j*beta)+phi))))
    F_2 = (1/k*R_2)*np.exp(-1j*np.cos((beta*C_2)/R_2)+1j*k*(Z_2 - C_2 * X_2)*(delta + 1j*beta))*(np.exp(1j*X_3*(k*C_2*(delta+1j*beta)+phi))-(np.exp(1j*X_2*(k*C_2*(delta+1j*beta)+phi))))
    F_3 = (1/k*R_3)*np.exp(-1j*np.cos((beta*C_3)/R_3)+1j*k*(Z_3 - C_3 * X_3)*(delta + 1j*beta))*(np.exp(1j*X_4*(k*C_3*(delta+1j*beta)+phi))-(np.exp(1j*X_3*(k*C_3*(delta+1j*beta)+phi))))
   
    print("F_0: {}".format(F_0))
    print("F_1: {}".format(F_1))
    print("F_2: {}".format(F_2))
    print("F_3: {}".format(F_3))
   
    
    n = (abs(F_0 + F_1 + F_2 + F_3)**2)*((d*np.cos(theta))**(-2))
    
    print("Diffraction Efficiency [n]: {}".format(n))
    
    return n

def gratingEfficiencyHARV(m,b,d,G = 0):
    ''''
    Return the theoretical diffraction efficiency of X-rays incident on an arbitrary grating
    with rectangular wire cross sections (method from Harvey et al (2019)):
        m: order of diffracted beam
        b: slit size
        d: grating periodicity
        G: grating type - 0 = amplitude grating, 1 = phase grating
    '''
    if G == 0:
        n = ((b**2)/(d**2))*(np.sinc((m*b)/d)**2)
        
    elif G == 1:
        if b != d/2:
            print("Error: Ratio of slit size to periodicity must be 1:2 for phase grating")
            import sys
            sys.exit()
        else :
            n = np.sinc(m/2)**2
    # print("absolute grating efficiency: {}".format(n))
    
    return n

def plotEfficiency():

    ampEff=[]
    phaEff=[]
    
    pEff_0 = []
    pEff_1 = []
    pEff_2 = []
    pEff_3 = []
    pEff_4 = []
    
    sEff_0 = []
    sEff_1 = []
    sEff_2 = []
    sEff_3 = []
    sEff_4 = []
    
    Eff_0 = []
    Eff_1 = []
    Eff_2 = []
    Eff_3 = []
    Eff_4 = []
    
    oEff = []
    
    
    period = np.array(range(1,200))
    # print("period: {}".format(period))
    slit = np.array(range(1,200))
    order = range(1,5)
    
    for p in period:
        E1 = gratingEfficiencyHARV(1,100e-9,1e-9*p,G=0)
        pEff_1.append(E1)
        E2 = gratingEfficiencyHARV(2,100e-9,1e-9*p,G=0)
        pEff_2.append(E2)
        E3 = gratingEfficiencyHARV(3,100e-9,1e-9*p,G=0)
        pEff_3.append(E3)
        E4 = gratingEfficiencyHARV(4,100e-9,1e-9*p,G=0)
        pEff_4.append(E4)
    # plt.plot(period, pEff_0, label = "m=0")
    plt.plot(period, pEff_1, label = "m=1")
    plt.plot(period, pEff_2, label = "m=2")
    plt.plot(period, pEff_3, label = "m=3")
    plt.plot(period, pEff_4, label = "m=4")
    plt.title("\u03B7\u2098 vs W\u209A (W\u209B=100 nm)")
    plt.xlabel("Period Width (W\u209A) [nm]")
    plt.ylabel("Grating Efficiency (\u03B7\u2098)")
    plt.legend()
    plt.show()
    
    
    for s in slit:
        # E0 = gratingEfficiencyHARV(0,1e-9*s,100e-9,G=0)
        # sEff_0.append(E0)
        E1 = gratingEfficiencyHARV(1,1e-9*s,100e-9,G=0)
        sEff_1.append(E1)
        E2 = gratingEfficiencyHARV(2,1e-9*s,100e-9,G=0)
        sEff_2.append(E2)
        E3 = gratingEfficiencyHARV(3,1e-9*s,100e-9,G=0)
        sEff_3.append(E3)
        E4 = gratingEfficiencyHARV(4,1e-9*s,100e-9,G=0)
        sEff_4.append(E4)
    # plt.plot(slit,sEff_0, label="m=0")
    plt.plot(slit,sEff_1, label="m=1")
    plt.plot(slit,sEff_2, label="m=2")
    plt.plot(slit,sEff_3, label="m=3")
    plt.plot(slit,sEff_4, label="m=4")
    plt.title("\u03B7\u2098 vs W\u209B (W\u209A=100 nm)")
    plt.xlabel("Slit Width (W\u209B) [nm]")
    plt.ylabel("Grating Efficiency (\u03B7\u2098)")
    plt.legend()
    plt.show()
    
    
    for s in slit:
        E0 = gratingEfficiencyHARV(0,1e-9*s,1e-9*p,G=0)
        Eff_0.append(E0)
        E1 = gratingEfficiencyHARV(1,1e-9*s,1e-9*p,G=0)
        Eff_1.append(E1)
        E2 = gratingEfficiencyHARV(2,1e-9*s,1e-9*p,G=0)
        Eff_2.append(E2)
        E3 = gratingEfficiencyHARV(3,1e-9*s,1e-9*p,G=0)
        Eff_3.append(E3)
        E4 = gratingEfficiencyHARV(4,1e-9*s,1e-9*p,G=0)
        Eff_4.append(E4)
        
    plt.plot(slit*0.005,Eff_0, label="m=0")
    plt.plot(slit*0.005,Eff_1, label="m=1")
    plt.plot(slit*0.005,Eff_2, label="m=2")
    plt.plot(slit*0.005,Eff_3, label="m=3")
    plt.plot(slit*0.005,Eff_4, label="m=4")
    plt.title("\u03B7\u2098 vs W\u209B/W\u209A")
    plt.yscale("log")
    plt.ylim(1e-4,1)
    plt.xlabel("W\u209B/W\u209A")
    plt.ylabel("Grating Efficiency (\u03B7\u2098)")
    plt.legend()
    plt.show()
    
    print("Shape of Eff_0: {}".format(np.shape(Eff_0)))
    print("Shape of Eff_1: {}".format(np.shape(Eff_1)))
    print("Shape of Eff_2: {}".format(np.shape(Eff_2)))
    print("Shape of Eff_3: {}".format(np.shape(Eff_3)))
    print("Shape of Eff_4: {}".format(np.shape(Eff_4)))
    
    tEff = np.array(Eff_0)+np.array(Eff_1)+np.array(Eff_2)+np.array(Eff_3)+np.array(Eff_4)
    totEff = np.array(Eff_1)+np.array(Eff_2)+np.array(Eff_3)+np.array(Eff_4)

    print("Shape of totEff: {}".format(np.shape(totEff)))
    
    plt.plot(slit*0.005,totEff, label = "m>0")
    plt.plot(slit*0.005,tEff, label = "m≥0")
    plt.plot(slit*0.005,Eff_0, label = "m=0")
    plt.title("\u03B7\u209C\u2092\u209C vs W\u209B/W\u209A")
    plt.xlabel("W\u209B/W\u209A")
    plt.ylabel("Total Grating Efficiency (\u03B7\u209C\u2092\u209C)")
    plt.legend()
    plt.show()
    
    
    
    for o in order:
        E = gratingEfficiencyHARV(o,1e-9,2e-9,G=1)
        oEff.append(E)
    ooEff=np.array([oEff[0], 
                   oEff[0]+oEff[1], 
                   oEff[0]+oEff[1]+oEff[2],  
                   oEff[0]+oEff[1]+oEff[2]+oEff[3]])
    print("Shape of oEff: {}".format(np.shape(oEff)))
    print("Shape of ooEff: {}".format(np.shape(ooEff)))
    plt.plot(order,oEff, label = "\u03B7\u2098")
    plt.plot(order,ooEff, label = "\u03B7\u209C\u2092\u209C")
    plt.title("\u03B7 vs m (phase grating)")
    plt.xlabel("Diffraction Order (m)")
    plt.ylabel("Grating Efficiency (\u03B7)")
    plt.legend()
    plt.show()
    
    
    
    # plt.plot(period, pEff, label="period (slit=100nm)")
    # plt.plot(slit, sEff, label="slit (period=200nm)")
    # plt.title("Grating Efficiency")
    # plt.xlabel("Size of slit/period [nm]")
    # plt.ylabel("Grating Efficiency (m=1)")
    # plt.legend()
    # plt.show()
    
    # for p in period:
    #     for s in slit:
    #         E = gratingEfficiencyHARV(1,1e-9*s,1e-9*p,G=0)
    #         ampEff.append(E)
            
    # print("Shape of ampEff: {}".format(np.shape(ampEff)))
    
    #     # s = np.reshape(stk.arS,(4,stk.mesh.nx,stk.mesh.ny))
    # plt.plot(ampEff)
    # plt.show()
    # # Eff = np.reshape(ampEff,(3,)))
    
    # # lamda = 6.7e-9
    # # p = 200e-9
    # # # m = 10
    # # NA = 0.1 #m*lamda/p
    # # # print("NA: {}".format(NA))
    
    # # m = NA*p/lamda
    # # print("m: {}".format(m))


def getImageData(filename):
    im = cv2.imread(filename, cv2.IMREAD_ANYDEPTH )   # open any bit depth image (more readable than -1)
    # show each image (can comment out for speed)
    #plt.imshow(im)
    #plt.show()
    
    # convert each image to array
    im = np.array(im)

    return im


def get_rect(x, y, width, height, angle):
    rect = np.array([(0, 0), (width, 0), (width, height), (0, height), (0, 0)])
    theta = (np.pi / 180.0) * angle
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    offset = np.array([x, y])
    transformed_rect = np.dot(rect, R) + offset
    
    return transformed_rect


def test():
    
    """ Testing theoretical efficiency """
    m=1                       #order of diffracted beam
    beta = 1.76493861*1e-3    #imaginary part of refractive index
    delta = 1-(2.068231*1e-2) #(1-real part) of refractive index
    d = 200e-9 #grating periodicity
    k = 197.3 * 6.7 #photon energy in units of (hbar*c)
    b = 72e-9
    theta = 0 #np.pi/2 #incident angle
    
    # print(" ")
    # print("-----Grating efficiency (DEL)-----")
    # gratingEfficiencyDEL(m, beta, delta, d, k, b,theta)
    
    print(" ")
    print("-----Grating efficiency (amplitude)-----")
    gratingEfficiencyHARV(1,100e-9,200e-9,G=0)
    
    print(" ")
    print("-----Grating efficiency (phase)-----")
    gratingEfficiencyHARV(1,100e-9,200e-9,G=1)


    """ Testing simulated efficiency """
    eMin = 1e8
    Nx = 150
    Ny = 150
    Nz = 1
    xMin = -10e-6
    xMax = 10e-6
    yMin = -10e-6
    yMax = 10e-6
    zMin = 100
    # Fx = 1/2
    # Fy = 1/2
    
    print(" ")
    print('Running Test:')
    print('building wavefront...')
    w = build_gauss_wavefront(Nx,Ny,Nz,eMin/1000,xMin,xMax,yMin,yMax,1,1e-6,1e-6,1) 
    
    wf0 = Wavefront(srwl_wavefront=w)
    
    """ Intensity from test Gaussian """
    # I = wf0.get_intensity()
    
    
    """ Intensity from tif file """
    I0 = getImageData('/home/jerome/WPG/intensity_atmask.tif')
    I1 = getImageData('/home/jerome/WPG/intensity_maskprop.tif') #getImageData('/home/jerome/WPG/intensityTot_maskprop.tif')
    
    I0_tot = np.sum(I0)
    I1_tot = np.sum(I1)
    
    Ir = I0_tot/I1_tot
    
    print("Intensity Ratio: {}".format(Ir))
    
    
    plt.imshow(I0)
    plt.title("at mask")
    plt.colorbar()
    plt.show()
    
    plt.imshow(I1)
    plt.title("after propagation")
    plt.colorbar()
    plt.show()
    
    Is = np.shape(I1)
    
    print("Shape of I: {}".format(Is))
    
    print(" ")
    print("-----Total Intensity-----")
    print("Itot: {}".format(I0_tot))

    """ Define region of interest to inspect separate orders """    
    ROI = ((int((Is[0]/2)-150),int(600)),
           (int((Is[0]/2)+150),int(1000)))
    
    x0,y0 = ROI[0][0], ROI[0][1]
    x1,y1 = ROI[1][0], ROI[1][1]
    
    
    A = I1[y0:y1,x0:x1]
    
    plt.imshow(A)
    plt.title('m=1')
    plt.colorbar()
    plt.show()
    
    Im_1 = np.sum(A)
    
    print("----- Intensity of m=+1-----")
    print("Im_1: {}".format(Im_1))
    
    # # Ri = (int((Is[0]/2)-150),int(600))
    # # Rf = (int((Is[0]/2)+150),int(900))
    
    # Ri = (int(1250),int(750))
    # Rf = (int(1250),int(750))

    # color = 100#(255,255,0)
    # thickness = 1000

    # Im = cv2.rectangle(I, 
    #               Ri, 
    #               Rf, 
    #               color, 
    #               thickness)
    
    # cv2.imshow('image', Im)
    
    # cv2.line(I,Ri,Rf,(0,250,0),100)


    # Convert the numpy array to an Image object.
    img = Image.fromarray(I1)

    # Draw a rotated rectangle on the image.
    draw = ImageDraw.Draw(img)
    rect = get_rect(x=750, y=750, width=500, height=500, angle=0.0)
    draw.polygon([tuple(p) for p in rect], fill='white')#, outline='white')
    # Convert the Image data to a numpy array.
    new_data = np.asarray(img)

    # Display the result using matplotlib.  (`img.show()` could also be used.)
    plt.imshow(new_data, cmap=plt.cm.gray)
    plt.show()


    E1 = Im_1/I1_tot
    
    print("Efficiency of m=1 order: {}".format(E1))
    
    
if __name__ == '__main__':
   test()