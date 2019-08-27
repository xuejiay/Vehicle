#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:41:24 2019

@author: xuejiaoyang
"""
import numpy as np
import matplotlib.pyplot as plt  
num_state=3
num_state_advance=4
num_state_advancez=10



step = 11*10
name = ['$x_e$ (cm)','$y_e$ (cm)','$\\theta_e$']
states=np.genfromtxt('States.csv', delimiter=',')
SDIbd=np.genfromtxt('StatesBD_SDI.csv', delimiter=',')
DIcst=np.genfromtxt('StatesBD_DI_CST.csv', delimiter=',')
DIcstz=np.genfromtxt('StatesBD_DI_CST_advance.csv', delimiter=',')


#for i in range(3):
#    F=plt.figure(i) 
#    for s in range(int(len(states)/step)):
#        sample = states[step*s:step*(s+1),:]
#        plt.plot(sample[:,0],sample[:,i+1],color = '0.75')
#    plt.plot(SDIbd[0:6,0],SDIbd[0:6,i+1],color = 'b')
#    plt.plot(SDIbd[0:6,0],SDIbd[0:6,i+num_state+1],color = 'b')
#    plt.plot(DIcst[:,0],DIcst[:,i+1],color = 'g')
#    plt.plot(DIcst[:,0],DIcst[:,i+num_state_advance+1],color = 'g')
#    plt.plot(DIcstz[:,0],DIcstz[:,i+1],color = 'm')
#    plt.plot(DIcstz[:,0],DIcstz[:,i+num_state_advancez+1],color = 'm')
#    plt.xlabel('s (m)',fontsize=16)
#    plt.ylabel(name[i],fontsize=16) 
#    plt.grid()
#    plt.tight_layout()
##    F.savefig(name[i],dpi=300)
    
    
sample_orig=np.genfromtxt('sample_orig.csv', delimiter=',')    
DI_CST_orig=np.genfromtxt('DI_CST_orig.csv', delimiter=',')
t = DIcst[:,0];
for i in range(3):
    F=plt.figure(i) 
    for s in range(int(len(states)/step)):
        plt.plot(t, sample_orig[:,s*3+i],color = '0.75')
    plt.plot(t,DI_CST_orig[:,i],color = 'g')
    plt.plot(t,DI_CST_orig[:,i+3],color = 'g')
    plt.plot(DIcstz[:,0],DIcstz[:,i+4+1],color = 'm')
    plt.plot(DIcstz[:,0],DIcstz[:,i+num_state_advancez+4+1],color = 'm')
    plt.xlabel('t (s)',fontsize=16)
    plt.ylabel(name[i],fontsize=16) 
    plt.grid()
    plt.tight_layout()
        
        
        
        
        
        