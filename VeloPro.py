# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 14:48:28 2022

@author: queco
"""
import math as m
import numpy as np
import os
from datetime import timedelta, date
import scipy
from scipy import stats
import matplotlib.pyplot as plt

def DCDFT(TimeInSec,RealVel,path1,saveFlag,Specname):
    if len(TimeInSec) <10:
        print('Too few points to perform DCDFT')
        Iw=[]
        omega=[]
    else:
        N = len(RealVel)
        mean = (sum(RealVel))/N
        for i in range(len(RealVel)):
            RealVel[i] -= mean
        if N % 2 == 0:
            offset = 0
        else:
            offset = 1
        wS = N/(max(TimeInSec)-min(TimeInSec))
        wN = wS/2
        dw = wS/N
        w = []
        W = -wN+(dw/2*offset)-dw
        while W <= wN-(dw/2*offset):
            W = W+dw
            w.append(W)
        win = []
        for n in range(0,len(RealVel)): 
            A = (46/(25*len(RealVel)))
            Window = A*((25/46)-((21/46)*m.cos((2*(m.pi)*n)/(len(RealVel)-1))))   
            win.append(Window)
        DataWithWin = [a*b for a,b in zip(RealVel,win)]
        a0 = N**(-1/2)
        ho = a0
        Fw = []
        #IWw2 = []
        for i in range(0,len(w)):
            w0 = w[i]
            FW = 0
            #tt = 0
            if w0 > 0.0:
                H0 = []
                H1 = []
                H2 = []
                h0 = []
                h1 = []
                h2 = []
                h1ih2 = []
                for j in range(0,len(TimeInSec)):
                    H_0 = 1
                    H_1 = m.cos(2*np.pi*w0*TimeInSec[j])
                    H_2 = m.sin(2*np.pi*w0*TimeInSec[j])
                    h_0 = a0*H_0
                    H0.append(H_0)
                    H1.append(H_1)
                    H2.append(H_2)
                    h0.append(h_0)
                IP_H0H1 = np.inner(H0,H1)
                IP_H0H2 = np.inner(H0,H2)
                IP_H1H1 = np.inner(H1,H1)
                IP_H1H2 = np.inner(H1,H2)
                IP_H2H2 = np.inner(H2,H2)
                IP_h0H1 = np.inner(h0,H1)
                IP_h0H2 = np.inner(h0,H2)
                a1 = ((IP_H1H1 - ((a0*IP_H0H1)**2))**(-(1/2)))
                a2 = ((IP_H2H2 - ((a0*IP_H0H2)**2) - ((a1*IP_H1H2)**2)
                    - ((a1*(a0**2)*IP_H0H1*IP_H0H2)**2)
                    + (2*((a1*a0)**2)*IP_H0H1*IP_H0H2*IP_H1H2))**(-(1/2)))
                for j in range(0,len(TimeInSec)):
                    h_1 = (a1*H1[j] - (a1*ho*IP_h0H1))
                    h1.append(h_1)
                IP_h1H2 = np.inner(h1,H2)
                for j in range(0,len(H1)):
                    h_2 = ((a2*H2[j]) - (a2*ho*IP_h0H2) - (a2*h1[j]*IP_h1H2))
                    h2.append(h_2)
                for j in range(0,len(TimeInSec)):
                    h1_ih2 = complex(h1[j],h2[j])
                    h1ih2.append(h1_ih2)
                    
                #c1 = [a*b for a,b in zip(DataWithWin,H1)]
                #c1=sum(c1)
                #c1=(a0*m.sqrt(2))*c1
                #print(c1)
                #c2 = [a*b for a,b in zip(DataWithWin,H2)]
                #c2=sum(c2)
                #c2=(a0*m.sqrt(2))*c2
                #tt=(c1**2)+(c2**2)
                #print(tt)
                RealVelComp = [a*b for a,b in zip(DataWithWin,h1ih2)]
                IP_DataWinh = sum(RealVelComp)
                FW = (IP_DataWinh)/(a0*m.sqrt(2))          
            Fw.append(FW)
            #IWw2.append(tt)
        #print(IWw2)    
        I = []
        for i in range(0,len(w)):
            II = 2*(a0**2)*(Fw[i]*np.conjugate(Fw[i]))
            I.append(II)
        wlog = []
        Ilog = []
        #I2log = []   #####
        for i in range(0,len(w)):
            Ww = w[i]    
            if Ww > 0:
                WW =m.log10(Ww)
                Ie = m.log10(I[i])
                #Ie2=m.log10(IWw2[i]) ####
                wlog.append(WW)
                Ilog.append(Ie)
                #I2log.append(Ie2) #####
        omega = np.array(wlog)
        Iw = np.array(Ilog)
        #Iw2 = np.array(I2log)
        if saveFlag == 1:
            SpecData=[(a,b) for a,b in zip(Iw,omega)]
            path2=(path1+"/Spectrum")
            os.chdir(path2)
            np.savetxt(Specname,SpecData)
        
    return Iw, omega#, Iw2

def loadVelo(path1,site,alt,yearstart,yearend,monthstart,monthend,D,
             saveFlag):
    os.chdir(path1)
    def datelist():
        start = date(int(yearstart), int(monthstart), int(1))
        end = date (int(yearend), int(monthend), int(1))
        for n in range(int((end-start).days)):
            yield (start + timedelta(n)).strftime('%Y%m%d')
    dates = list(datelist())
    DataInd = []
    VelRaw = []
    for i in range(0,len(dates)):
        file = [str(dates[int(i)]) +'.' + site + '.' + alt
            + D + '_1min.txt']
        if os.path.isfile(str(file[0])) is True:
            dataind = (np.loadtxt(str(file[0]),delimiter=',',usecols=(9,)))
            velraw = (np.loadtxt(str(file[0]),delimiter=',',usecols=(8,)))
            VelRaw = np.concatenate((VelRaw, velraw), axis=0)
            DataInd = np.concatenate((DataInd, dataind), axis=0)
        else:
            fakedataind = np.zeros((1440))
            fakevelraw = np.zeros((1440))
            VelRaw = np.concatenate((VelRaw, fakevelraw), axis=0)
            DataInd = np.concatenate((DataInd, fakedataind), axis=0)
    mnt = -1
    TimeInSec = []
    RealVel = []
    for veltrue, t in zip(VelRaw,DataInd): 
        if t == 1:
            mnt = mnt+1
            sec = mnt*60
            TimeInSec.append(sec)
            RealVel.append(veltrue)
        if t == 0:
            mnt = mnt+1
    if not RealVel:
        RealVel=[1,3,2,3,1,3]
        TimeInSec=[10,1000,4000,15000,60000,10000]
    if len(RealVel)<=4:
        RealVel=[1,3,2,3,1,3]
        TimeInSec=[10,1000,4000,15000,60000,10000]
    if saveFlag == 1:
        VelData=[(a,b) for a,b in zip(RealVel,TimeInSec)]
        path2=(path1+"/Velocities")
        os.chdir(path2)
        np.savetxt(site+str(yearstart)+monthstart+alt+D,VelData)
        
    return RealVel, TimeInSec

def loadSpec(path1,Specname):
    path2=(path1+"/Spectrum")
    os.chdir(path2)
    if os.path.isfile(Specname) is True:
        Iw1 = np.array((np.loadtxt((Specname),delimiter=' ',usecols=(0,))))
        Omega1 = np.array((np.loadtxt((Specname),delimiter=' ',usecols=(1,))))
    else:
        Iw1 = [1,2]
        Omega1 = [1,3]
    if len(Iw1)<10:
        print('Empty Spectrum for: '+Specname)
        Iw = []
        Omega = []
    else:
        Iw = []
        Omega = []
        for i in range(0,len(Iw1)):
            if m.isnan(Iw1[i]) is False:
                Iw.append(Iw1[i])
                Omega.append(Omega1[i])
    return Iw, Omega

def func(Omega,a,b,c,d):
    return a/(1+(np.e**((Omega-b)/c)))+d

def fitSpec(Iw,Omega,cutOff,plotTog):
    if len(Iw)<10:
        print('too few points to fit slope')
        fitResults = (np.nan,np.nan,np.nan,np.nan,np.nan)
    else:
        try:
            popt, pcov = scipy.optimize.curve_fit(func,Omega,Iw,
                            bounds = ([1,-6,.1,-10],[4,-4,.5,-4]))
        except:
            popt=[0.999999,1,1.0001]
            print('Warning: FD fit did not work') #has never happened
            

        cutOff2 = 1-cutOff
        x1 = popt[2]*np.log((1/cutOff2)-1)+popt[1]
        #if x1<=-5.2:
        #    x1=-5.2
        x2 = popt[2]*np.log((1/cutOff)-1)+popt[1]
        
        if plotTog ==1:
            plt.figure(1)
            plt.plot(Omega,Iw,'k',Omega, func(Omega,*popt),'r')
            plt.xlim([-6,-3.5])
            plt.figure(1)
            plt.plot(x1,func(x1,*popt),'b*',x2,func(x2,*popt),'b*',markersize=20)
        
        Islope = []
        Oslope = []
        for i in range(0,len(Omega)):
            O = Omega[i]
            P = Iw[i]
            if O >= x1 and O <= x2:
                Oslope.append(O)
                Islope.append(P)
        if len(Islope)<3:
            print('did not use, Too few points for slope')
            fitResults = (np.nan,np.nan,np.nan,np.nan,np.nan)
        else:
            slope, intercept, r_value, p_value, std_err = stats.linregress(Oslope,Islope)
            fitResults = (slope , std_err , r_value , x1 , x2)
            
            if plotTog ==1:
                X=np.arange(x1-1,x2+1,.05)
                Y=slope*X+intercept
                
                plt.figure(1)
                plt.plot(X,Y,'g^')
                plt.xlim([-6,-3.25])
                plt.grid('on',which='major', axis='both')
                plt.show(1)
            
    return fitResults
    