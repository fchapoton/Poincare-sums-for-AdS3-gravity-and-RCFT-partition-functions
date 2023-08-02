#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:14:29 2020

@author: viraj29
"""

import numpy as np
from sage import *


#-----------------Functions Required----------------------#

def hdim(lamb):
    res1 = (np.power(lamb,2) - 1)/(4*(k_lev+2))
    return res1

def conj(vec):
    res2 = np.zeros((len(vec),len(vec)),dtype=complex)
    for i in range(0,len(vec)):
        for j in range(0,len(vec)):
            res2[i,j] = (np.conj(vec[i]))*(vec[j])
    return res2

def tns(n):
    tnmat = np.zeros((dim,dim),dtype=complex)
    smat = np.zeros((dim,dim), dtype=complex)
    for lamb in range(1,dim+1):
        tnmat[lamb-1,lamb-1] = np.exp( np.complex(0,2*pi*n*hdim(lamb)) )
    for lamb in range(1,dim+1):
        for rho in range(1,dim+1):
            smat[lamb-1,rho-1] = (np.sqrt(2/(k_lev + 2)))*np.sin((pi*lamb*rho)/(k_lev+2))
    res3 = np.matmul(tnmat,smat)
    return res3

def tmat_b(n):
    return np.array([[1,n],[0,1]])

smat_b = np.array([[0,-1],[1,0]])

def word_list(trial,res4):
    if trial[1,0] == 0:
        sgn = trial[0,0]
        res4 = np.append(res4, [sgn*trial[0,1], sgn])
        return res4
    else:
        if abs(trial[0,0])>=abs(trial[1,0]):
            rem = np.mod(trial[0,0],trial[1,0])
            div = (trial[0,0] - rem)/trial[1,0]
            res4 = np.append(res4,-div)
            return word_list(np.matmul(tmat_b(-div),trial),res4)
        else:
            if len(res4)==1:
                res4 = [1]
            else:
                res4 = np.append(res4,1)
            return word_list(np.matmul(smat_b,trial),res4)

def word_fin(vec1):
    if vec1[0]==1:
        res5 = [1]
    else:
        res5 = [0]
    for i in range(1,len(vec1)-2,2):
        res5 = np.append(res5, [-vec1[i]])
    res5 = np.append(res5, [vec1[-2]])
    return res5

def poincare_sum(mat):
    for lamb in range(1,dim+1):
        for rho in range(1,dim+1):
            if mat[lamb-1,rho-1] != 0:
                file.write(str(mat[lamb-1,rho-1]) + "chi(" + str(lamb) +  ") * chi(" +  str(rho) +  ") \n")

#-------------------------------Code----------------------------------#

file = open("WZW_SU(2)_check","w+")

for k_lev in range(1,6):
    file.write("****************************************************************************************************************************************************** \n")
    file.write("****************************************************************************************************************************************************** \n")
    pi=np.pi
    # k_lev = int (input("Enter SU(2) level = ")) # SU(2)_k
    dim = k_lev + 1 
    n0 = 4*(k_lev + 2) # semi-conductor
    gc = Gamma1(n0)
    coset_gc = FareySymbol(gc).coset_reps()
    coset_gc_n = np.zeros((len(coset_gc),2,2))
    
    for i in range(0,len(coset_gc)):
        for j in range(0,2):
            for k in range(0,2):
                coset_gc_n[i,j,k] = coset_gc[i][j,k]
    
    
    smatrix = tns(0)
    tns_mat = np.zeros((n0,dim,dim),dtype=complex)
    for i in range(0,n0):
        tns_mat[i,:,:]=tns(i)
    
    file.write("SU(2) Level = " + str(k_lev) + " | Semi-conductor N0 = " + str(n0) + " | Central Charge = " + str((3*k_lev)/(k_lev+2)) + "\n")
    file.write(str(gc)+ "\n")
    file.write("---------------------------------- \n")
    
    print("SU(2) Level = " + str(k_lev) + " | Semi-conductor N0 = " + str(n0) + " | Central Charge = " + str((3*k_lev)/(k_lev+2)) + "\n")
    print(str(gc)+ "\n")
    print("---------------------------------- \n")
    
    
    nta=np.zeros((len(coset_gc_n),dim,dim), dtype='complex')
    
    for i in range(0,len(coset_gc_n)):
        word = word_fin(word_list(coset_gc_n[i],[0]))
        if word[0]==1:
            fmat = tns_mat[0,:,:]
        else:
            fmat = np.identity(dim)
        for j in range(1,len(word)-1):
            fmat = np.matmul(fmat, tns_mat[int(word[j]) , : , : ])
        fmat = np.matmul(fmat, np.matmul(tns_mat[int(word[-1]), : , : ], tns_mat[ 0 , : , :]))
        nta[i,:,:] = fmat
    
    
    #--------------------------------Output-----------------------------------#
    fres=np.zeros((dim,dim,dim), dtype='complex')
    for lamb in range(1,dim+1):
        for i in range(0,len(coset_gc_n)):
            fres[lamb-1,:,:] = fres[lamb-1,:,:] + conj(np.transpose(nta[i,:,lamb-1]))
    
    fres1=np.round(fres,decimals=3)
    
    file.write("Poincare sum: P(1) \n")
    poincare_sum(fres1[0])
    
    for rho in range(2,dim+1):
        itr=1
        for lamb in range(1,rho):
            if np.array_equal(fres1[lamb-1], fres1[rho-1]):
                file.write("\n ")
                file.write("P(" + str(rho) + ") = P(" + str(lamb) + ") \n")
                itr=0
                break
        if itr==1:
            file.write("---------------------------------- \n")
            file.write("Poincare sum: P(" + str(rho) + ") \n")
            poincare_sum(fres1[rho-1])

file.write("\n\n")
file.close()
