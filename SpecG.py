# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:29:56 2017

@author: meisu
"""

import numpy as np 
from sys import argv 
from dump import dump 
from scipy.fftpack import fft, fftfreq
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

script, input_file = argv 

if rank == 0:
    
    correlation_length = int(raw_input("Enter the length of correlation:"))
    cross_section = float(raw_input("Enter the cross-section of interface (A^2):"))
    dt = float(raw_input("Enter the timestep used in MD (in ps):"))
    
    #Reading the dump files
    d = dump(input_file)
    time_steps = d.time()

    atom_list = d.vecs(time_steps[0],"id")
    atom_list = np.array(atom_list, dtype = np.int16) 
    num_atoms = len(atom_list)

    del_T = dt*(time_steps[1] - time_steps[0])
    scaling_const = 2/(cross_section*del_T*num_atoms)
    
    
else: 
    correlation_length = None 
    cross_section = None 
    dt = None 
    num_atoms = None 
    scaling_const = None
    atom_list = None 
    d = None 
#%%
#Sharing the input settings among processors  
correlation_length = comm.bcast(correlation_length, root = 0)
cross_section = comm.bcast(cross_section, root = 0)
dt = comm.bcast(dt, root = 0)
num_atoms = comm.bcast(num_atoms, root = 0)
scaling_const = comm.bcast(scaling_const, root = 0) 
d = comm.bcast(d, root = 0)
#%%
#Allocating memory
q_12x = np.zeros(correlation_length, dtype = np.complex64)
q_12y = np.zeros(correlation_length, dtype = np.complex64)
q_12z = np.zeros(correlation_length, dtype = np.complex64)
#%%
atom_list_recvbuf = np.empty(num_atoms/size, dtype = np.int16)

#%%
comm.Scatter(atom_list, atom_list_recvbuf, root = 0)
atom_list = atom_list_recvbuf
#%%
# Debugging 
for i in range(size): 
    if i == rank: 
        print "Atom list - proc:" + str(i) + " \n" 
        print atom_list

#%%
#Obtaining velocities, performing fourier transform
for i in atom_list:
    
    if i%1000 ==0:
        print "Number of atoms covered:"+str(i) + "(proc:" + str(rank)+")"
    
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

    
        q_12x[j] += np.sum(scaling_const * FX[:n_sample] *
                            np.conj(VX[j:n_sample+j]))/n_sample
        q_12y[j] += np.sum(scaling_const * FY[:n_sample]
                            * np.conj(VY[j:n_sample+j]))/n_sample
        q_12z[j] += np.sum(scaling_const * FZ[:n_sample] 
                            * np.conj(VZ[j:n_sample+j]))/n_sample

Q_tot_12 = q_12x + q_12y + q_12z

comm.Barrier()
#%%                
Q_tot_12_global = np.zeros(len(Q_tot_12), dtype = np.complex64)

comm.Reduce(Q_tot_12, [Q_tot_12_global,MPI.COMPLEX], 
		op = MPI.SUM, root = 0)
  
np.savetxt('Q_tot_12.data',Q_tot_12_global)
 
        
        
    
