# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:29:56 2017

@author: meisu
"""

import numpy as np 
from sys import argv 
from dump import dump 
from scipy.fftpack import fft, fftfreq

script, input_file1, input_file2 = argv 

correlation_length = int(raw_input("Enter the length of correlation:"))
cross_section = float(raw_input("Enter the cross-section of interface (A^2):"))
dt = float(raw_input("Enter the timestep used in MD (in ps):"))
#%%
#Reading the dump files
d1 = dump(input_file1)
d2 = dump(input_file2)

time_steps1 = d1.time()
time_steps2 = d2.time()

num_atoms1 = len(d1.vecs(time_steps1[0],"vx"))
num_atoms2 = len(d2.vecs(time_steps2[0],"vx"))

del_T1 = dt*(time_steps1[1] - time_steps1[0])
del_T2 = dt*(time_steps2[1] - time_steps2[0]) 
#%%
print   "Collecting velocities and force data..."

#Allocating memeory for velocities      
vx1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
vy1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
vz1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)  


#Allocating memory for forces
fx1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
fy1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
fz1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)

#Obtaining velocities and force data 
for i in range(len(time_steps1)):
    
    vx1[:,i], vy1[:,i], vz1[:,i], fx1[:,i], fy1[:,i], fz1[:,i] = d1.vecs(
                                            time_steps1[i], "vx", "vy", "vz", 
                                            "c_F_C_M[1]", "c_F_C_M[2]", "c_F_C_M[3]")

del d1 
#%%
FX1 = []
FY1 = []
FZ1 = []

VX1 = [] 
VY1 = []
VZ1 = []

for i in range(num_atoms1):
    FX1.append(fft(fx1[i,:]))
    FY1.append(fft(fy1[i,:]))
    FZ1.append(fft(fz1[i,:]))

    VX1.append(fft(vx1[i,:]))
    VY1.append(fft(vy1[i,:]))
    VZ1.append(fft(vz1[i,:]))
    
FX1 = np.array(FX1)
FY1 = np.array(FY1)
FZ1 = np.array(FZ1)

VX1 = np.array(VX1)
VY1 = np.array(VY1)
VZ1 = np.array(VZ1)

del vx1,vy1,vz1,fx1,fy1,fz1

#%%
#Allocating memeory for velocities              
vx2 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)   
vy2 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)
vz2 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)


#Allocating memory for forces
fx2 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)
fy2 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)
fz2 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)




for i in range(len(time_steps2)): 
    
    vx2[:,i], vy2[:,i], vz2[:,i], fx2[:,i], fy2[:,i], fz2[:,i] = d2.vecs(
                                            time_steps2[i], "vx", "vy", "vz", 
                                            "c_F_M_C[1]", "c_F_M_C[2]", "c_F_M_C[3]")
del d2
#%%
FX2 = []
FY2 = []
FZ2 = []

VX2 = [] 
VY2 = []
VZ2 = []

for i in range(num_atoms2):
    FX2.append(fft(fx2[i,:]))
    FY2.append(fft(fy2[i,:]))
    FZ2.append(fft(fz2[i,:]))

    VX2.append(fft(vx2[i,:]))
    VY2.append(fft(vy2[i,:]))
    VZ2.append(fft(vz2[i,:]))
    
FX2 = np.array(FX2)
FY2 = np.array(FY2)
FZ2 = np.array(FZ2)

VX2 = np.array(VX2)
VY2 = np.array(VY2)
VZ2 = np.array(VZ2)

del vx2,vy2,vz2,fx2,fy2,fz2    
#%%
print "Performing cross-correlation....."

q_12x = np.array(correlation_length)
q_12y = np.array(correlation_length)
q_12z = np.array(correlation_length)

scaling_const = 2/(cross_section*del_T1)

for j in range(correlation_length): 
    
    n_sample = len(FX1)-j

    for i in range(n_sample): 
        
        
        q_12x[j] = (scaling_const * FX1[j] * np.conj(VX1[i+j]))/n_sample
        q_12y[j] = (scaling_const * FY1[j] * np.conj(VY1[i+j]))/n_sample
        q_12z[j] = (scaling_const * FZ1[j] * np.conj(VZ1[i+j]))/n_sample

Q_tot_12 = q_12x + q_12y + q_12z
                
q_21x = np.array(correlation_length)
q_21y = np.array(correlation_length)
q_21z = np.array(correlation_length)

scaling_const = 2/(cross_section*del_T2)

for j in range(correlation_length): 
    
    n_sample = len(FX2)-j

    for i in range(n_sample): 
        
        
        q_21x[j] = (scaling_const * FX2[j] * np.conj(VX2[i+j]))/n_sample
        q_21y[j] = (scaling_const * FY2[j] * np.conj(VY2[i+j]))/n_sample
        q_21z[j] = (scaling_const * FZ2[j] * np.conj(VZ2[i+j]))/n_sample

Q_tot_21 = q_21x + q_21y + q_21z
        
        
    