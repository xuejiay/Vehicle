#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:13:58 2019

@author: xuejiaoyang
"""


from interval import interval, imath
import numpy as np


    
def erroPos2origPos(pe,pref,flag)  :
    
    Y0=pe[0]
    Y1=pe[1]
    Y2=pe[2]
    Y3=pe[3]
    Y4=pe[4]
    phid=pref[2]
    x1d=pref[0]
    x2d=pref[1]
    l=2.0
    if flag == 0:
        p_orig=np.zeros(5)
        p_orig[0] = np.cos(phid)*Y0-np.sin(phid)*Y1+x1d
        p_orig[1] = np.sin(phid)*Y0+np.cos(phid)*Y1+x2d
        p_orig[2] = Y2+phid
        p_orig[4] = Y3
        p_orig[3] = np.arctan(l*Y4)
    else:
        p_orig=[0.]*5
        p_orig[0] = imath.cos(phid)*Y0-imath.sin(phid)*Y1+x1d
        p_orig[1] = imath.sin(phid)*Y0+imath.cos(phid)*Y1+x2d
        p_orig[2] = Y2+phid
        p_orig[4] = Y3
        p_orig[3] = imath.atan(l*Y4)  
        
    return p_orig 




num_state_constr=10

#xs = np.genfromtxt('ref.csv', delimiter=',')
states=np.genfromtxt('States.csv', delimiter=',')
SDI=np.genfromtxt('StatesBD_SDI.csv', delimiter=',')
DI=np.genfromtxt('StatesBD_DI_CST.csv', delimiter=',')
DI_advanced=np.genfromtxt('StatesBD_DI_CST_advance.csv', delimiter=',')
t=DI_advanced[:,0]
step = 8*10


ref = DI_advanced[:,14+1:17+1]


# original coordinates of sampled states
state_origcoor = 5
sample_orig = np.zeros([len(t),64*state_origcoor])
for s in range(64):          
    sample = states[step*s:step*(s+1),1:6] 
    for i in range(len(t)):
        sample_orig[i,s*state_origcoor:(s+1)*state_origcoor] = erroPos2origPos(sample[i],ref[i],0)


# original coorinates of DI with constraints
DI_CST_orig = np.zeros([len(t),state_origcoor*2])
for i in range(len(t)):
    DI_e = [0.]*state_origcoor
    ref_interval = [0.]*state_origcoor
    for j in range(state_origcoor):
        DI_e[j] = interval([DI[i,j+1],DI[i,j+num_state_constr+1]])
    for k in range(3):
        ref_interval[k] = interval(ref[i,k])
        
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
