# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 14:13:33 2026

@author: cohen
"""

import matplotlib.pyplot as plt 
import numpy as np
import scipy
from tqdm import tqdm

plt.rcParams['font.size'] = 16


wL=0.014 #is 3.25 micron
E0=-0.0046#-0.008#-0.0046
cycle=2*np.pi/wL
t_max=cycle/5
cmap = plt.colormaps['viridis'] 



Eg=0.11760#3.2 eV in aau
a=5.32
kappa=2*np.pi/a
coeff_v = [-0.0928,0.0705,0.0200,-0.0012,0.0029,0.0006]
coeff_c = [0.0898,-0.0814,-0.0024,-0.0048,-0.0003,-0.0009]

Ev = lambda k: sum(coeff_v[i]*np.cos(i*k*a) for i in range(6))
Ec = lambda k: Eg+sum(coeff_c[i]*np.cos(i*k*a) for i in range(6))

Ev_diff = lambda k: sum(-i*a*coeff_v[i]*np.sin(i*k*a) for i in range(6))
Ec_diff = lambda k: sum(-i*a*coeff_c[i]*np.sin(i*k*a) for i in range(6))

laser = lambda t: E0*np.cos(wL*t)

energies = []
t0s=[]
trs=[]


plt.figure()
plt.title("x for t")
plt.xlabel("t in laser cylces")
plt.ylabel("x in grid units")
for t0 in tqdm(np.linspace(-t_max, t_max, 100)):
    k= lambda t: -E0/wL*(np.sin(wL*t)-np.sin(wL*t0)) #sol to k_dot=-E(t) , k(t0)=0
    
    delta_t=0.1
    
    k_t=[0]
    xv=[0]
    xc=[0]
    x_diff=[0]
    t=[t0]
    epsilon=0.1
    t_lim=1.5*cycle
    start=True
    while (start or abs(xc[-1]-xv[-1])>epsilon) and t[-1]<t_lim:
        if start and abs(xc[-1]-xv[-1])>epsilon:
            start=False
            
        t.append(t[-1]+delta_t)
        k_cur=k(t[-1])
        k_t.append(k_cur)
        xv.append(xv[-1]+delta_t*Ev_diff(k_cur))
        xc.append(xc[-1]+delta_t*Ec_diff(k_cur))
        # t.append(t[-1]+delta_t)
        x_diff.append(xc[-1]-xv[-1])
    
    if t[-1]<t_lim:
        energies.append(Ec(k(t[-1]))-Ev(k(t[-1])))
        t0s.append(t0)
        trs.append(t[-1])
    color=cmap((t0+t_max)/(2*t_max))
    
    plt.plot(np.array(t)/cycle,np.array(x_diff)/a,color=color)
    # plt.plot(np.array(t)/cycle,k_t)#,color=color)
    
plt.show()


# plt.figure()
# plt.title("energy for t")


# t_E = np.linspace(-t_max,max(trs),1000)

# plt.plot(t_E,laser(t_E),label="E")
# plt.plot(t0s,energies,'.',label="t0")
# plt.plot(trs,energies,'.',label="t return")
# #plt.plot(t_diffs,energies,'.',label="t return - t0")
# plt.legend()
# plt.show()
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