#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:14:29 2020

@author: viraj29
"""

import numpy as np
from sage import *


#-----------------Functions Required----------------------#

def hdim(a1):
    res1 = (a1*(N1-a1))/(2*N1)
    return res1

def conj(vec): # Changing done

    #--------------Reduced space--------------------
    
    res2 = np.zeros((dim_new,dim_new),dtype=complex)
    for i in range(0,len(vec)):
        if i>dim_new-1:
            i1=N1-i
        else:
            i1=i
        for j in range(0,len(vec)):
            if j>dim_new-1:
                j1=N1-j
            else:
                j1=j
            res2[i1,j1] = res2[i1,j1] + (np.conj(vec[i]))*(vec[j])        
    return res2

    #--------------Full space--------------------
    
    # res2 = np.zeros((len(vec),len(vec)),dtype=complex)
    # for i in range(0,len(vec)):
    #     for j in range(0,len(vec)):
    #         res2[i,j] = (np.conj(vec[i]))*(vec[j])
    # return res2
    
def tns(n):
    tnmat = np.zeros((dim,dim),dtype=complex)
    smat = np.zeros((dim,dim), dtype=complex)
    for a1 in range(0,dim):
        tnmat[a1,a1] = np.exp( np.complex(0,2*pi*n*hdim(a1)) )  # w/o the -c/24 term
    for a1 in range(0,dim):
        for b1 in range(0,dim):
            smat[a1,b1] = (1/np.sqrt(N1))*np.exp( np.complex(0,(-2*pi*a1*b1)/N1) )
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

def poincare_sum(mat): # Changing now
    for a1 in range(0,dim_new):
        for b1 in range(0,dim_new):
            if mat[a1,b1] != 0:
                file.write(str(mat[a1,b1]) + "chi(" + str(a1) +  ") * chi(" +  str(b1) +  ") \n")

#-------------------------------Code----------------------------------#

file = open("WZW_SU(N=1_to_10)_1.txt","w+")

for N1 in range(1,11):
    file.write("****************************************************************************************************************************************************** \n")
    file.write("****************************************************************************************************************************************************** \n")
    pi=np.pi
    dim = N1 
    if np.mod(N1,2) != 0:
        dim_new = int((N1+1)/2)
    if np.mod(N1,2) == 0:
        dim_new = int((N1+2)/2)
    n0 = int(2*N1) # semi-conductor
    gc = Gamma1(n0)
    coset_gc = FareySymbol(gc).coset_reps()
    coset_gc_n = np.zeros((len(coset_gc),2,2))
    
    for i in range(0,len(coset_gc)):
        for j in range(0,2):
            for k in range(0,2):
                coset_gc_n[i,j,k] = coset_gc[i][j,k]
    
    smatrix = tns(0)
    smatrix3 = np.matmul( np.matmul(tns(0),tns(0)) , tns(0))
    tns_mat = np.zeros((n0,dim,dim),dtype=complex)
    for i in range(0,n0):
        tns_mat[i,:,:]=tns(i)
    
    file.write("SU(" + str(N1) + ")1 | Semi-conductor N0 = " + str(n0) + " | Central Charge = " + str(N1-1) + "\n")
    file.write(str(gc)+ "\n")
    file.write("---------------------------------- \n")
    
    print("SU(N = " + str(N1) + ")1 | Semi-conductor N0 = " + str(n0) + " | Central Charge = " + str(N1-1) + "\n")
    print(str(gc)+ "\n")
    print("---------------------------------- \n")
    
    
    nta=np.zeros((len(coset_gc_n),dim,dim), dtype='complex')
    
    for i in range(0,len(coset_gc_n)):
        word_sgn = word_list(coset_gc_n[i],[0])[-1]
        word = word_fin(word_list(coset_gc_n[i],[0]))
        if word[0]==1:
            fmat = -tns_mat[0,:,:] # change to -1 done
        else:
            fmat = np.identity(dim)
        for j in range(1,len(word)-1):
            fmat = -np.matmul(fmat, tns_mat[int(word[j]) , : , : ]) # change to -1 done
        fmat = word_sgn*np.matmul(fmat, np.matmul(tns_mat[int(word[-1]), : , : ], smatrix3 )) # word_sgn done needed  # tns_mat[ 0 , : , :]
        nta[i,:,:] = fmat
    
    # for i in range(0,len(coset_gc_n)):
    #     word = word_fin(word_list(coset_gc_n[i],[0]))
    #     if word[0]==1:
    #         fmat = tns_mat[0,:,:]
    #     else:
    #         fmat = np.identity(dim)
    #     for j in range(1,len(word)-1):
    #         fmat = np.matmul(fmat, tns_mat[int(word[j]) , : , : ])
    #     fmat = np.matmul(fmat, np.matmul(tns_mat[int(word[-1]), : , : ], tns_mat[ 0 , : , :]))
    #     nta[i,:,:] = fmat
    
    
    #--------------------------------Output-----------------------------------#
    fres=np.zeros((dim,dim_new,dim_new), dtype='complex') # Changing done
    for a1 in range(0,dim):
        for i in range(0,len(coset_gc_n)):
            fres[a1,:,:] = fres[a1,:,:] + conj(np.transpose(nta[i,:,a1]))
    
    fres1=np.round(fres,decimals=3)
    
    file.write("Poincare sum: P(0) \n")
    poincare_sum(fres1[0])
    
    for a1 in range(1,dim):
        itr=1
        for b1 in range(0,a1):
            if np.array_equal(fres1[a1], fres1[b1]):
                file.write("\n ")
                file.write("P(" + str(a1) + ") = P(" + str(b1) + ") \n")
                itr=0
                break
        if itr==1:
            file.write("---------------------------------- \n")
            file.write("Poincare sum: P(" + str(a1) + ") \n")
            poincare_sum(fres1[a1])

file.write("\n\n")
file.close()
