#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:13:58 2019

@author: xuejiaoyang
"""


import numpy as np

#import Dubin_full as case
from scipy.integrate import odeint
#from coordtrans2 import rhs
import matplotlib.pyplot as plt   


num_state_advance=2
xs = np.genfromtxt('ref.csv', delimiter=',')

SDI = np.genfromtxt('origPosSDI.csv', delimiter=',')
DI = np.genfromtxt('origPosDI_CST.csv', delimiter=',')

SDI_s = np.genfromtxt('SDI_s.csv', delimiter=',')
DI_s = np.genfromtxt('DIconstraints_s.csv', delimiter=',')
advDI_s = np.genfromtxt('advanceDI_s.csv', delimiter=',')


xs = np.genfromtxt('ref.csv', delimiter=',')
X_meth3 = np.genfromtxt('StatesBD_SDI.csv', delimiter=',')
el =X_meth3[:,1]
eu =X_meth3[:,1+num_state_advance]
t=X_meth3[:,0]


origPos = np.zeros([len(X_meth3),4])
   
for i in range(len(X_meth3)):
    theta_orth = xs[i,1]+np.pi/2
    x = xs[i,2]
    y = xs[i,3]
    xel = el[i]*np.cos(theta_orth)
    yel = el[i]*np.sin(theta_orth)
    xl = x+xel
    yl = y+yel
    
    xeu = eu[i]*np.cos(theta_orth)
    yeu = eu[i]*np.sin(theta_orth)
    xu = x+xeu
    yu = y+yeu  
    
    origPos[i,0] = xl
    origPos[i,1] = yl
    origPos[i,2] = xu
    origPos[i,3] = yu
    
    

step = 801*2    
states = np.genfromtxt('States.csv', delimiter=',')  
origPos = np.zeros([step,2])   
F =plt.figure(0)
for s in range(int(len(states)/step)):           
    sample = states[step*s:step*(s+1),:]  
    el =sample[:,1]
    for i in range(len(sample)):
        theta_orth = xs[i,1]+np.pi/2
        x = xs[i,2]
        y = xs[i,3]
        xel = el[i]*np.cos(theta_orth)
        yel = el[i]*np.sin(theta_orth)
        xl = x+xel
        yl = y+yel

        
        origPos[i,0] = xl
        origPos[i,1] = yl

    plt.plot(origPos[:,0],origPos[:,1],color = '0.75')
plt.plot(DI[:,0],DI[:,1],color = 'g')    
plt.plot(DI[:,0+2],DI[:,1+2],color = 'g')  
plt.plot(SDI[:,0],SDI[:,1],color = 'b')    
plt.plot(SDI[:,0+2],SDI[:,1+2],color = 'b')  

plt.plot(SDI_s[0:7,4],SDI_s[0:7,5],color = 'r')  
plt.plot(SDI_s[0:7,4+6],SDI_s[0:7,5+6],color = 'r')  
plt.plot(DI_s[0:7,4],DI_s[0:7,5],color = 'm')  
plt.plot(DI_s[0:7,4+7],DI_s[0:7,5+7],color = 'm')  
plt.plot(advDI_s[0:7,4],advDI_s[0:7,5],color = 'k')  
plt.plot(advDI_s[0:7,4+9],advDI_s[0:7,5+9],color = 'k')  



#plot for shor horizon
#plt.plot(DI[0:50,0],DI[0:50,1],color = 'g')    
#plt.plot(DI[0:50,0+2],DI[0:50,1+2],color = 'g')  
#plt.plot(SDI[:,0],SDI[:,1],color = 'b')    
#plt.plot(SDI[:,0+2],SDI[:,1+2],color = 'b')  
#
#plt.plot(SDI_s[0:7,4],SDI_s[0:7,5],color = 'r')  
#plt.plot(SDI_s[0:7,4+6],SDI_s[0:7,5+6],color = 'r')  
#plt.plot(DI_s[0:7,4],DI_s[0:7,5],color = 'm')  
#plt.plot(DI_s[0:7,4+7],DI_s[0:7,5+7],color = 'm')  
#plt.plot(advDI_s[0:7,4],advDI_s[0:7,5],color = 'k')  
#plt.plot(advDI_s[0:7,4+9],advDI_s[0:7,5+9],color = 'k')  
#


name = ['x (m)','y (m)']
plt.xlabel(name[0],fontsize=16)
plt.ylabel(name[1],fontsize=16)  
F.savefig('path',dpi=300)
    
    
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
