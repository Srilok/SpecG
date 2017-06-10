# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:29:56 2017

@author: meisu
"""

import numpy as np 
from sys import argv 
from dump import dump 

script, input_file1, input_file2 = argv 

#%%
#Reading the dump files
d1 = dump(input_file1)
d2 = dump(input_file2)

time_steps1 = d1.time()
time_steps2 = d2.time()

num_atoms1 = len(d1.vecs(time_steps1[0],"fx"))
num_atoms2 = len(d2.vecs(time_steps2[0],"fx"))

#%%
print   "Collecting velocities and force data..."

#Allocating memeory for velocities      
vx1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
vy1 = np.zeros([num_atoms2,len(time_steps1)], dtype = np.float32)
vz1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)  


#Allocating memory for forces
fx1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
fy1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
fz1 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)

#Obtaining velocities and force data 
for i in range(len(time_steps1)):
    
    vx1[:,i], vy1[:,i], vz1[:,i], fx1[:,i], fy1[:,i], fz1[:,i] = d1.vecs(
                                            time_steps1[i], "vx", "vy", "vz", 
                                            "fx", "fy", "fz")

del d1 
#%%

#Allocating memeory for velocities              
vx2 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)   
vy2 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
vz2 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)


#Allocating memory for forces
fx2 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
fy2 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)
fz2 = np.zeros([num_atoms1,len(time_steps1)], dtype = np.float32)




for i in range(len(time_steps2)): 
    
    vx2[:,i], vy2[:,i], vz2[:,i], fx2[:,i], fy2[:,i], fz2[:,i] = d2.vecs(
                                            time_steps2[i], "vx", "vy", "vz", 
                                            "fx", "fy", "fz")
del d2

#%%

#Computing F_i from Q_1->2

Fi_12 = np.zeros([num_atoms1], dtype = np.float32) 

print    "Computing F_i (= sum(f_ij, j belongs to material 2))..."

for i in range(num_atoms1): 
    for j in range(num_atom2): 
        Fi_12[i] += 
        

print    "Performing fourier transform..."


    