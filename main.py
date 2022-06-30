# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 13:56:10 2022

@author: queco
"""
#The function of this main code is to allow the user to specify which radar
#they want to process and any height range etc. After that the code will take
#over.
import os
import VeloPro
import matplotlib.pyplot as plt
import numpy as np
SITE = ['aumond','egbert','eureka','gananoque','harrow','markstay','mcgill'
      ,'negrocreek','walsingham','wilberforce','abitibi'];
mon = ['01','02','03','04','05','06','07','08','09','10','11','12']
RDir = ['N','E','S','W']
MonthStart = ['01','02','03','04','05','06','07','08','09','10','11','12']
MonthEnd   = ['02','03','04','05','06','07','08','09','10','11','12','01']
Alt2 = ['14.0','13.5','13.0','12.5','12.0','11.5','11.0','10.5','10.0','9.5'
                    ,'9.0','8.5','8.0','7.5','7.0','6.5','6.0','5.5','5.0',
                    '4.5','4.0','3.5','3.0','2.5','2.0','1.5','1.0']
pwd = os.getcwd()
saveSpec = 0 #keep 1 if you want to save the Velocity series
saveVel = 0 #keep 1 if you want to save the Spectrum series
DirOfData='C:/Users/queco/Desktop/'
for Site in range(1,2):
    site = SITE[Site]
    for year in range(2011,2012):
        YearStart  = [year,year,year,year,year,year,year,year,year,year,year,
                      year]
        YearEnd  = [year,year,year,year,year,year,year,year,year,year,year,
                    (year+1)]
        for month in range(0,5):
            yearstart  = YearStart[month]
            yearend    = YearEnd[month]
            monthstart = MonthStart[month]
            monthend   = MonthEnd[month]
            for Dir in range(2,3):
                D = RDir[Dir] #Direction
                path1 = (DirOfData+site+"Unpack"+"/re"+str(year)+monthstart+"."+
                         site+"/"+D)
                try:
                    os.chdir(path1)
                    for ALT1 in range(10,11):
                        alt = Alt2[ALT1]
                        print(site,str(year),monthstart,RDir[Dir],alt)
                        [RealVel,TimeInSec] = VeloPro.loadVelo(path1,site,
                                 alt,yearstart,yearend,monthstart,monthend,D,
                                 saveVel)
                        print(len(TimeInSec),'number of real time data points')
                        plt.plot(TimeInSec,RealVel,'b.')
                        plt.show()
                        
                        Specname=(site+str(year)+monthstart+alt+D)

                        #Iw , omega = VeloPro.DCDFT(TimeInSec,RealVel,path1,saveSpec,Specname)
                        #plt.plot(omega,Iw)
                        #plt.xlim((-7,-3.5))
                        #plt.show()
                        
                        Iw , Omega = VeloPro.loadSpec(path1,Specname)
                        #plt.plot(Omega,Iw)
                        #plt.xlim((-7,-3.5))
                        #plt.show()
                        
                        cutOff=0.15
                        plotTog=1
                        fitResults = VeloPro.fitSpec(Iw,Omega,cutOff,plotTog)
                        
                        
                except:
                    pass
os.chdir(pwd)