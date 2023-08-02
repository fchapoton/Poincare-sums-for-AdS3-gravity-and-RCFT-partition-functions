#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:14:29 2020

@author: viraj29
"""

import numpy as np

#-------------Changes to be made before start------------#
#                 Values of (p,q) (lines 18-19)
#                 lcm of denominators of h_{r,s} (lines 20)
#                 Seed primary    (lines 92) 
#--------------------------------------------------------#

pi=np.pi
p=11 # Minimal model
q=12
hlcm=528
dim=int((p-1)*(q-1)/2) # dimension of rep.
nta = np.zeros((1,dim),dtype=complex)
dnta=1
ll=0
rdim = q-2 # (Check) Needs to be changed for different models


#---------------------------Functions--------------------------------#

def hdim(r,s):
    res1 = (np.power((p*r) - (q*s),2) - np.power(p-q,2) )/(4*p*q)
    return res1
    
def conj(vec):
    res2 = np.zeros((len(vec),len(vec)),dtype=complex)
    for i in range(0,len(vec)):
        for j in range(0,len(vec)):
            res2[i,j] = (np.conj(vec[i]))*(vec[j])
    return res2

def tns(n):
    tnmat = np.zeros((dim,dim),dtype=complex)
    for r in range(1,rdim+1):
        for s in range(r,rdim+1):
            for c in range(1,rdim+1):
                for d in range(c,rdim+1):
                    i1 = ((r-1)*rdim) - (r*(r-1)/2) + s - 1
                    i2 = ((c-1)*rdim) - (c*(c-1)/2) + d - 1
                    if i1==i2:
                        tnmat[int(i1),int(i2)] = np.exp(complex(0,2*pi*n*hdim(r,s)))

    smat = np.zeros((dim,dim), dtype=complex)
    
    for r in range(1,rdim+1):            
        for s in range(r,rdim+1):
            for c in range(1,rdim+1):
                for d in range(c,rdim+1):
                    i1 = ((r-1)*rdim) - (r*(r-1)/2) + s - 1
                    i2 = ((c-1)*rdim) - (c*(c-1)/2) + d - 1
                    smat[int(i1),int(i2)] = np.sqrt(8/(p*q))*(np.power(-1,1+(s*c)+(r*d)))*(np.sin(pi*p*r*c/q))*(np.sin(pi*q*s*d/p))
    
    res3 = np.matmul(tnmat,smat)
    return res3

def phase(vec1,vec2,nitr):
    if np.abs(vec1[nitr])>1e-6:
        if np.abs(vec2[nitr])<1e-6:
            return 0
        else:
            return np.angle(vec2[nitr]) - np.angle(vec1[nitr])
    else:
        return phase(vec1,vec2,nitr+1)

def new_contr(trial1,mitr,nta1):
    if mitr < len(nta1):
        check = np.transpose(nta1[mitr,:])
        phi = phase(trial1,check,0)
        if np.all(np.abs(trial1 - np.multiply(np.exp(complex(0,-phi)),check)) < 1e-6):
            return 0
        else:
            return new_contr(trial1,mitr+1,nta1)
    else:
        return 1


smatrix = np.round(tns(0),decimals=3)
tns_mat = np.zeros((hlcm,dim,dim),dtype=complex)
for i in range(0,hlcm):
    tns_mat[i,:,:]=tns(i)

#-----------------Non-trivial operations-------------------#
nta[0,8] = 1.0 + 0j # seed primary (Range - 0 to dim-1)
print("Initial Seed: ",nta)
new=1
while (new != 0):
    ul=dnta
    new=0
    for n in range(0,hlcm): # depeds upon lcm of h_{r,s} - Value of N_{0}
        for tn in range(ll,ul):
            trial = np.matmul(tns_mat[n,:,:],np.transpose(nta[tn,:]))
            itr=1
            for mitr in range(0,dnta):
                check = np.transpose(nta[mitr,:])
                phi = phase(trial,check,0)
                if np.all(np.abs(trial - np.multiply(np.exp(complex(0,-phi)),check)) < 1e-6):
                    itr=0
                    break
            if itr==1:
                nta = np.append(nta, [trial], axis=0)
                dnta = dnta + 1
                new = new + 1
    ll=ul
    print("New =", new)

#---------Final Result----------#
fres=0
fres1=0
for l in range(0,len(nta[:,0])):
    fres = fres + conj(nta[l,:])
    
fres1 = np.round(fres,decimals=3)
print(fres1)

for i in range(0,dim):
    for j in range(0,dim):
        if fres1[i,j] != 0:
            print(fres1[i,j],"\chi[", i, "] * \chi[", j, "]")
            


#----------------------------------#
#itr=new_contr(trial,0,nta)
#
#
#            itr=1
#            for mitr in range(0,dnta):
#                check = np.transpose(nta[mitr,:])
#                phi = phase(trial,check,0)
#                if np.all(np.abs(trial - np.multiply(np.exp(complex(0,-phi)),check)) < 1e-6):
#                    itr=0





