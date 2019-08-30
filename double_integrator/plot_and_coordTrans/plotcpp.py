#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:41:24 2019

@author: xuejiaoyang
"""
import numpy as np
import matplotlib.pyplot as plt  
num_state=7
num_state_advance=10
num_state_advancez=17


#step = 8*10
#name = ['$e_t$ (m)','$e_n$ (m)', '$e_{\phi}$','v', '$\kappa_{\delta}$']
#states=np.genfromtxt('States.csv', delimiter=',')
#SDIbd=np.genfromtxt('StatesBD_SDI.csv', delimiter=',')
#DIcst=np.genfromtxt('StatesBD_DI_CST.csv', delimiter=',')
#DIcstz=np.genfromtxt('StatesBD_DI_CST_advance.csv', delimiter=',')
#ref=np.genfromtxt('ref.csv', delimiter=',')
#for i in range(num_state):
#    F=plt.figure(i) 
##    plt.plot(ref[:,0],ref[:,i+2],color = 'b')
#    for s in range(64):
##    for s in range(int(len(states)/step)):
#        sample = states[step*s:step*(s+1),:]
##        for i in range(2):
#        plt.plot(sample[:,0],sample[:,i+1],color = '0.75')
#    plt.plot(SDIbd[:,0],SDIbd[:,i+1],color = 'b')
#    plt.plot(SDIbd[:,0],SDIbd[:,i+num_state+1],color = 'b')
#    plt.plot(DIcst[:,0],DIcst[:,i+1],color = 'g')
#    plt.plot(DIcst[:,0],DIcst[:,i+num_state_advance+1],color = 'g')
#    plt.plot(DIcstz[:,0],DIcstz[:,i+1],color = 'm')
#    plt.plot(DIcstz[:,0],DIcstz[:,i+num_state_advancez+1],color = 'm')
#    plt.xlabel('t (s)',fontsize=16)
#    plt.ylabel(name[i],fontsize=16) 
#    plt.grid()
#    plt.tight_layout()
#    F.savefig(name[i],dpi=300)
    
    
DIcstz=np.genfromtxt('StatesBD_DI_CST_advance.csv', delimiter=',')    
sample_orig=np.genfromtxt('sample_orig.csv', delimiter=',')    
DI_CST_orig=np.genfromtxt('DI_CST_orig.csv', delimiter=',')
name = ['x (cm)','y (cm)','$\\phi$', '$\\delta$']
t = DIcstz[:,0];
for i in range(5):
    F=plt.figure(i) 
    for s in range(64):
        plt.plot(t, sample_orig[:,s*5+i],color = '0.75')
    plt.plot(t,DI_CST_orig[:,i],color = 'g')
    plt.plot(t,DI_CST_orig[:,i+5],color = 'g')
    plt.plot(DIcstz[:,0],DIcstz[:,i+10+1],color = 'm')
    plt.plot(DIcstz[:,0],DIcstz[:,i+num_state_advancez+10+1],color = 'm')
    plt.xlabel('t (s)',fontsize=16)
    plt.ylabel(name[i],fontsize=16) 
    plt.grid()
    plt.tight_layout()