# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 16:42:50 2026

@author: cohen
"""

import matplotlib.pyplot as plt 
import numpy as np
import scipy
from tqdm import tqdm

plt.rcParams['font.size'] = 16


wL=0.00227 #is two micron
E0=-0.006
cycle=2*np.pi/wL
t_max=cycle/5
cmap = plt.colormaps['viridis'] 


plt.figure()
plt.title("X for Time")
plt.xlabel("Time in laser cylces")
plt.ylabel("X (a.u.)")
for t0 in tqdm(np.linspace(-t_max, t_max, 40)):
    t0+=t_max/2
    x = lambda t: E0/wL**2*(np.cos(wL*t)-np.cos(wL*t0))+E0/wL*np.sin(wL*t0)*(t-t0)
    
    
    t=np.linspace(t0,2*cycle,1000)

    for i,tt in enumerate(t[10:]):
        if x(tt)<0:
            tf=scipy.optimize.fsolve(x,tt)#t[:i+10+1]
            t=np.linspace(t0,tf,1000)
            break
    
    color=cmap((t0+t_max/2)/(3/2*t_max))
    plt.plot(np.array(t)/cycle,x(t),color=color)
    
plt.show()



energies = []
t0s=[]
trs=[]
t_diffs=[]
for t0 in tqdm(np.linspace(-t_max, t_max, 1000)):
    x = lambda t: E0/wL**2*(np.cos(wL*t)-np.cos(wL*t0))+E0/wL*np.sin(wL*t0)*(t-t0)
    v = lambda t: -E0/wL*np.sin(wL*t)+E0/wL*np.sin(wL*t0)
    
    cycle=2*np.pi/wL
    t=np.linspace(t0,2*cycle,1000)

    if x(t[1])>0:
        for tt in t:
            if x(tt)<0:
                tf=scipy.optimize.fsolve(x,tt)#t[:i+10+1]
                t0s.append(t0)
                trs.append(tf)
                t_diffs.append(tf-t0)
                energies.append(v(tf)**2/2)
                break
            
# plt.figure()
# plt.title("energy for t")

# E=lambda t: -E0/wL*np.cos(wL*t)
# t_E = np.linspace(-t_max,max(trs),1000)

# plt.plot(t_E,E(t_E),label="E")
# plt.plot(t0s,energies,'.',label="t0")
# plt.plot(trs,energies,'.',label="t return")
# #plt.plot(t_diffs,energies,'.',label="t return - t0")
# plt.legend()
# plt.show()
laser=lambda t: -E0/wL*np.cos(wL*t)
fig, ax1 = plt.subplots()
plt.title("Recombination Energy vs Time")

t_E = np.linspace(-t_max, max(trs), 1000)

color1 = 'tab:red'
ax1.set_xlabel("Time (a.u.)")
ax1.set_ylabel("Energy (a.u.)", color=color1)
ax1.plot(t0s, energies, '.', color='blue', label="t0 (Birth)",markersize=9)
ax1.plot(trs, energies, '.', color='red', label="t return",markersize=9)
ax1.tick_params(axis='y', labelcolor=color1)
ax1.legend(loc='upper left')

ax2 = ax1.twinx()  
color2 = 'tab:blue'
ax2.set_ylabel("Electric Field", color=color2)  
ax2.plot(t_E, laser(t_E), '-', color=color2, alpha=0.5, label="E-field")
ax2.tick_params(axis='y', labelcolor=color2)

plt.show()