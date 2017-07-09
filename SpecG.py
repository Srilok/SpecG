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
q_12x = np.zeros(correlation_length, dtype = np.complex64)
q_12y = np.zeros(correlation_length, dtype = np.complex64)
q_12z = np.zeros(correlation_length, dtype = np.complex64)
#%%
print   "Collecting velocities and force data..."

#Obtaining velocities, performing fourier transform
for i in atom_list:
    print i
    
    if i%1000 ==0:
        print "Number of atoms covered:"+str(i)
    
    vx, vy, vz, fx, fy, fz = d.atom(i, "vx", "vy", "vz",
                                   "c_F_C_M[1]", "c_F_C_M[2]", "c_F_C_M[3]") 
    
    VX = fft(vx)
    VY = fft(vy)
    VZ = fft(vz)
    
    FX = fft(fx)
    FY = fft(fy)
    FZ = fft(fz)


    del vx,vy,vz,fx,fy,fz
    #Performing cross-correlation for each atom and accumulating the sum    
    for j in range(correlation_length): 
    
        n_sample = len(FX)-j

    
        q_12x[j] += np.sum(scaling_const * FX[:n_sample] 
                            * np.conj(VX[j:n_sample+j]))/n_sample
        q_12y[j] += np.sum(scaling_const * FY[:n_sample]
                            * np.conj(VY[j:n_sample+j]))/n_sample
        q_12z[j] += np.sum(scaling_const * FZ[:n_sample] 
                            * np.conj(VZ[j:n_sample+j]))/n_sample

Q_tot_12 = q_12x + q_12y + q_12z
                
#        
        
    