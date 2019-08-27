#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 16:41:24 2019

@author: xuejiaoyang
"""
import numpy as np
import matplotlib.pyplot as plt  
num_state=2
num_state_advance=3
num_state_advancez=5
#df = pd.read_csv("States.csv")
#with open('States.csv', newline='') as csvfile:
#    data = list(csv.reader(csvfile))
#s = 78
#
#print(data[101*s:101*(s+1)])
step = 101*2
name = ['e (m)','$\\theta_e$']
#states=np.genfromtxt('States.csv', delimiter=',')
#SDIbd=np.genfromtxt('StatesBD_SDI2.csv', delimiter=',')
DIcst=np.genfromtxt('StatesBD_DI_CST.csv', delimiter=',')
DIcstz=np.genfromtxt('StatesBD_DI_CST_advance.csv', delimiter=',')
ref=np.genfromtxt('ref.csv', delimiter=',')
for i in range(2):
    F=plt.figure(i) 
#    plt.plot(ref[:,0],ref[:,i+2],color = 'b')
#    for s in range(int(len(states)/step)):
##    for s in range(30):
#        sample = states[step*s:step*(s+1),:]
#        for i in range(2):
#            plt.plot(sample[:,0],sample[:,i+1],color = '0.75')
#    plt.plot(SDIbd[0:30,0],SDIbd[0:30,i+1],color = 'b')
#    plt.plot(SDIbd[0:30,0],SDIbd[0:30,i+num_state+1],color = 'b')
    plt.plot(DIcst[:,0],DIcst[:,i+1],color = 'g')
    plt.plot(DIcst[:,0],DIcst[:,i+num_state_advance+1],color = 'g')
#    plt.plot(DIcstz[:,0],DIcstz[:,i+1],color = 'm')
#    plt.plot(DIcstz[:,0],DIcstz[:,i+num_state_advancez+1],color = 'm')
    plt.xlabel('s (m)',fontsize=16)
    plt.ylabel(name[i],fontsize=16) 
    plt.grid()
    plt.tight_layout()
#    F.savefig(name[i],dpi=300)