# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:29:56 2017

@author: meisu
"""

import numpy as np 
from sys import argv 
from dump import dump 
from scipy.fftpack import fft, fftfreq

script, input_file = argv 

correlation_length = int(raw_input("Enter the length of correlation:"))
cross_section = float(raw_input("Enter the cross-section of interface (A^2):"))
dt = float(raw_input("Enter the timestep used in MD (in ps):"))
#%%
#Reading the dump files
d = dump(input_file)
#d2 = dump(input_file2)

time_steps = d.time()
#%%
#time_steps2 = d2.time()
atom_list = d.vecs(time_steps[0],"id")
atom_list = np.array(atom_list, dtype = np.int16) 
num_atoms = len(atom_list)
#num_atoms2 = len(d2.vecs(time_steps2[0],"vx"))

del_T = dt*(time_steps[1] - time_steps[0])
#del_T2 = dt*(time_steps2[1] - time_steps2[0]) 

scaling_const = 2/(cross_section*del_T*num_atoms)
#%%
#Allocating memory
K_12x = np.zeros(correlation_length, dtype = np.complex64)
K_12y = np.zeros(correlation_length, dtype = np.complex64)
K_12z = np.zeros(correlation_length, dtype = np.complex64)

q_12x = np.zeros(correlation_length, dtype = np.complex64)
q_12y = np.zeros(correlation_length, dtype = np.complex64)
q_12z = np.zeros(correlation_length, dtype = np.complex64)
#%%
#Obtaining frequency list 

vx = d.atom(1,"vx")
freq_list = fftfreq(correlation_length,0.01*1e-12)
freq_positive = freq_list[0:len(freq_list)/2]
#%%
print   "Collecting velocities and force data..."

#Obtaining velocities, performing fourier transform
for i in atom_list:
    
    if i%1000 ==0:
        print "Number of atoms covered:"+str(i)
    
    
    vx, vy, vz, fx, fy, fz = d.atom(i, "vx", "vy", "vz",
                                   "c_F_C_M[1]", "c_F_C_M[2]", "c_F_C_M[3]")
    
    vx=np.array(vx)
    vy=np.array(vy)
    vz=np.array(vz)
    
    fx=np.array(fx)
    fy=np.array(fy)
    fz=np.array(fz)
    

    #Performing cross-correlation for each atom and accumulating the sum    
    for j in range(correlation_length): 
    
        n_sample = len(vx)-j

    
        K_12x[j] += np.sum(scaling_const * fx[:n_sample] 
                            * vx[j:n_sample+j])/n_sample
        K_12y[j] += np.sum(scaling_const * fy[:n_sample]
                            * vy[j:n_sample+j])/n_sample
        K_12z[j] += np.sum(scaling_const * fz[:n_sample] 
                            * vz[j:n_sample+j])/n_sample

    del vx,vy,vz,fx,fy,fz 
    
    q_12x += 2*np.real(fft(K_12x))
    q_12y += 2*np.real(fft(K_12y))
    q_12z += 2*np.real(fft(K_12z))


Q_tot_12 = q_12x + q_12y + q_12z
#%%

        
    
