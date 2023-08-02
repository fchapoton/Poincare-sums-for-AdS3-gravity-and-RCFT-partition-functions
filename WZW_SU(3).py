# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 14:52:18 2020

@author: palash

Created on Thu Oct 15 13:14:29 2020

@author: viraj29
"""

from sage import *
import numpy as np

pi=np.pi

#-----------------Functions Required----------------------#

def hdim(lamb1,lamb2):
    res1 = (np.power(lamb1,2) + np.power(lamb2,2) + (lamb1*lamb2) - 3)/(3*(k_lev+3))
    return res1

def conj(vec): # NEW CHANGES MADE
    res2 = np.zeros((dim_new,dim_new),dtype=complex)
    for i in range(1,len(vec)+1):
        x1 = lu_name(i)
        if x1[0]>x1[1]:
            xi1 = x1[0] + (x1[1]-1)*(2+k_lev-x1[1]) #x1[0] + ( (k_lev+1)*(k_lev+2) - (k_lev+2-x1[1])*(k_lev+3-x1[1])  )/2
        else:
            xi1= x1[1] + (x1[0]-1)*(2+k_lev-x1[0])
        for j in range(1,len(vec)+1):
            y1 = lu_name(j)
            if y1[0]>y1[1]:
                yi1 = y1[0] + (y1[1]-1)*(2+k_lev-y1[1]) #y1[0] + ( (k_lev+1)*(k_lev+2) - (k_lev+2-y1[1])*(k_lev+3-y1[1])  )/2
            else:
                yi1 = y1[1] + (y1[0]-1)*(2+k_lev-y1[0])
            res2[int(xi1)-1,int(yi1)-1] = res2[int(xi1)-1,int(yi1)-1] + (np.conj(vec[i-1]))*(vec[j-1])
    return res2

def enx(j1,j2,j3,j4,lamb1,lamb2,mu1,mu2): # not have omitted the '-i' overall factor
    return (complex(0,-1)/(np.sqrt(3)*(k_lev+3)))*np.exp(  complex(0,((-2*pi)/(3*(k_lev+3)))*( j1*lamb1*mu1 + j2*lamb1*mu2 + j3*lamb2*mu1 + j4*lamb2*mu2 ) )  )

def tns(n):
    tnmat = np.zeros((dim,dim),dtype=complex)
    smat = np.zeros((dim,dim), dtype=complex)    
    i1=0
    for lamb1 in range(1,k_lev+2):
        for lamb2 in range(1,k_lev-lamb1+3): # CHANGE----------
            i2=0
            for mu1 in range(1,k_lev+2):
                for mu2 in range(1,k_lev-mu1+3): # CHANGE----------
                    if i1==i2:
                        tnmat[int(i1),int(i2)] = np.exp(complex(0,2*pi*n*hdim(lamb1,lamb2)))
                    i2=i2+1
            i1=i1+1
    i1=0
    for lamb1 in range(1,k_lev+2):
        for lamb2 in range(1,k_lev-lamb1+3): # CHANGE----------
            i2=0
            for mu1 in range(1,k_lev+2):
                for mu2 in range(1,k_lev-mu1+3): # CHANGE----------
                    smat[int(i1),int(i2)] = enx(2,1,1,2,lamb1,lamb2,mu1,mu2) + enx(-1,-2,1,-1,lamb1,lamb2,mu1,mu2) + enx(-1,1,-2,-1,lamb1,lamb2,mu1,mu2) - enx(-1,-2,-2,-1,lamb1,lamb2,mu1,mu2) - enx(2,1,1,-1,lamb1,lamb2,mu1,mu2) - enx(-1,1,1,2,lamb1,lamb2,mu1,mu2)
                    i2=i2+1
            i1=i1+1
                                       
    res3 = np.matmul(tnmat,smat)
    return res3

#                    i1 =  -np.power(lamb1,2) + ((k_lev+3)*lamb1) + lamb2-k_lev-3-1  # (lamb1 - 1)*(k_lev + 2 - (0.5*lamb1)) + lamb2 - 1 
#                    i2 =  -np.power(mu1,2) + ((k_lev+3)*mu1) + mu2-k_lev-3-1        # (mu1 - 1)*(k_lev + 2 - (0.5*mu1)) + mu2 - 1


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
        res5 = [1]  # may need to change to -1
    else:
        res5 = [0]
    for i in range(1,len(vec1)-2,2):
        res5 = np.append(res5, [-vec1[i]])  # maybe already has -1
    res5 = np.append(res5, [vec1[-2]])
    return res5
    
def lu_name(num):
    num_count=1
    for lamb1 in range(1,k_lev+2):
        for lamb2 in range(1,k_lev-lamb1+3):             # CHANGE
            if num_count==num:
                return [lamb1,lamb2]
            else:
                num_count=num_count+1

def lu_name_new(num):
    num_count=1
    for lamb1 in range(1,k_lev+2):
        for lamb2 in range(lamb1,k_lev-lamb1+3):        # CHANGE
            if num_count==num:
                return [lamb1,lamb2]
            else:
                num_count=num_count+1

def poincare_sum(mat):
    for i in range(0,dim_new):
        for j in range(0,dim_new):
            if mat[i,j] != 0:
                file.write(str(mat[i,j]) + "chi" + str(lu_name_new(i+1)) + " * chi" + str(lu_name_new(j+1)) + "\n")

#-------------------------------Code----------------------------------#

file = open("WZW_SU(3)_1_to_10.txt","w+")

for k_lev in range(1,11):
    file.write("************************************************************************************************************************************ \n")
    file.write("************************************************************************************************************************************ \n")
    dim = int((k_lev+1)*(k_lev+2)/2)
    if np.mod(k_lev,2) != 0:
        dim_new = int((k_lev+1)*(k_lev+3)/4)
    if np.mod(k_lev,2) == 0:
        dim_new = int((k_lev+2)*(k_lev+2)/4)
    n0 = 3*(k_lev + 3)                                 # semi-conductor
    gc = Gamma1(n0)                                     # no CHANGE
    coset_gc = FareySymbol(gc).coset_reps()
    coset_gc_n = np.zeros((len(coset_gc),2,2))
    
    for i in range(0,len(coset_gc)):
        for j in range(0,2):
            for k in range(0,2):
                coset_gc_n[i,j,k] = coset_gc[i][j,k]
        
    smatrix3 = np.matmul( np.matmul(tns(0),tns(0)) , tns(0))
    tns_mat = np.zeros((n0,dim,dim),dtype=complex)
    for i in range(0,n0):
        tns_mat[i,:,:]=tns(i)
    
    file.write("SU(3) Level = " + str(k_lev) + " | Dimension =" + str(dim_new) + " | Semi-conductor N0 = " + str(n0) + " | Central Charge = " + str((8*k_lev)/(k_lev+3)) + "\n")
    file.write(str(gc)+ "\n")
    file.write("---------------------------------- \n")
    
    print("SU(3) Level = " + str(k_lev) + " | Dimension = " + str(dim_new) + " | Semi-conductor N0 = " + str(n0) + " | Central Charge = " + str((8*k_lev)/(k_lev+3)) + "\n")
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
    
    
    #--------------------------------Output-----------------------------------#
    fres=np.zeros((dim,dim_new,dim_new), dtype='complex') # NEW CHANGES MADE
    for lamb in range(1,dim+1):
        for i in range(0,len(coset_gc_n)):
            fres[lamb-1,:,:] = fres[lamb-1,:,:] + conj(np.transpose(nta[i,:,lamb-1]))
    
    fres1=np.round(fres,decimals=3)
    
    file.write("Poincare sum: P(1,1) \n")
    poincare_sum(fres1[0])
    
    for rho in range(2,dim+1):
        itr=1
        for lamb in range(1,rho):
            if np.array_equal(fres1[lamb-1], fres1[rho-1]):
                file.write("\n ")
                file.write("P(" + str(lu_name(rho)) + ") = P(" + str(lu_name(lamb)) + ") \n")
                itr=0
                break
        if itr==1:
            file.write("---------------------------------- \n")
            file.write("Poincare sum: P(" + str(lu_name(rho)) + ") \n")
            poincare_sum(fres1[rho-1])

file.write("\n\n")
file.close()
