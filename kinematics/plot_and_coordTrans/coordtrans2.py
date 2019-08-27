#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:13:58 2019

@author: xuejiaoyang
"""


from interval import interval, imath
import numpy as np


    
def erroPos2origPos(pe,pref,flag)  :
    
    xr=pref[0]
    yr=pref[1]
    thetar=pref[2]
    xe=pe[0]
    ye=pe[1]
    thetae=pe[2]
    thetac=thetar-thetae
    
    
    if flag == 0:
        pc=np.zeros(3)
        pc[1] = yr-np.sin(thetac)*xe-np.cos(thetac)*ye
        pc[0] = xr-np.cos(thetac)*xe+np.sin(thetac)*ye
        pc[2] = thetac
    else:
        pc=[0]*3
        pc[1] = yr-imath.sin(thetac)*xe-imath.cos(thetac)*ye
        pc[0] = xr-imath.cos(thetac)*xe+imath.sin(thetac)*ye
        pc[2] = thetac    
    return pc




num_state_constr=4

#xs = np.genfromtxt('ref.csv', delimiter=',')
states=np.genfromtxt('States.csv', delimiter=',')
SDI=np.genfromtxt('StatesBD_SDI.csv', delimiter=',')
DI=np.genfromtxt('StatesBD_DI_CST.csv', delimiter=',')
DI_advanced=np.genfromtxt('StatesBD_DI_CST_advance.csv', delimiter=',')
t=DI_advanced[:,0]
step = 11*10
name = ['x (cm)','y (cm)','$\\theta$']

ref = DI_advanced[:,7+1:10+1]


# original coordinates of sampled states
state_origcoor = 3
sample_orig = np.zeros([len(t),int(len(states)/step)*state_origcoor])
for s in range(int(len(states)/step)):          
    sample = states[step*s:step*(s+1),1:4] 
    for i in range(len(t)):
        sample_orig[i,s*state_origcoor:(s+1)*state_origcoor] = erroPos2origPos(sample[i],ref[i],0)


# original coorinates of DI with constraints
DI_CST_orig = np.zeros([len(t),state_origcoor*2])
for i in range(len(t)):
    DI_e = [0.]*state_origcoor
    ref_interval = [0.]*state_origcoor
    for j in range(state_origcoor):
        DI_e[j] = interval([DI[i,j+1],DI[i,j+num_state_constr+1]])
        ref_interval[j] = interval(ref[i,j])
        
    x_orig = erroPos2origPos(DI_e,ref_interval,1)
    for l in range(state_origcoor):
        DI_CST_orig[i,l] = x_orig[l][0][0]
        DI_CST_orig[i,l+state_origcoor] = x_orig[l][0][1]
        
np.savetxt('DI_CST_orig.csv', DI_CST_orig, delimiter=',')    
np.savetxt('sample_orig.csv', sample_orig, delimiter=',')          
        
#for j in range(2):
#    F =plt.figure(j)
##    for s in range(num_sample):    
##        plt.plot(t, xs[s,:,j+2],color = '0.75')
##    plt.grid()
##    plt.tight_layout()
##    plt.xlabel('s (m)',fontsize=16)
##    plt.ylabel('x (m)',fontsize=16)      
#
##plt.plot(t, X_meth1[:,i],color = 'b')
##plt.plot(t, X_meth1[:,i+nx],color = 'b')
##    plt.plot(t, X_meth2[:,i],color = 'r')
##    plt.plot(t, X_meth2[:,i+case.paramAIV.nx],color = 'r') 
##    plt.plot(t, origPos[:,j],color = 'm')
##    plt.plot(t, origPos[:,j+3],color = 'm') 
#    plt.plot(t, origPos[:,j],color = 'g')
#    plt.plot(t, origPos[:,j+2],color = 'g')  
##    plt.plot(t, origPos2[:,j],color = 'm')
##    plt.plot(t, origPos2[:,j+2],color = 'm')  
#    plt.xlabel('s (m)',fontsize=16)
#    plt.ylabel(name[j],fontsize=16)  
##plt.plot(t, origPos[:,1],color = 'g')
##plt.plot(t, origPos[:,4],color = 'g')  
##    F.savefig(name[j],dpi=300)   
#
##np.savetxt('origPos.csv', origPos, delimiter=',') 
#
