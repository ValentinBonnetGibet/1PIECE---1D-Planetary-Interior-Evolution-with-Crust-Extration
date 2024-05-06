#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 09:49:07 2022

@author: valentinbonnetgibet
"""

import numpy as np
import math as m

nproc = 1024

N = int(m.log2(nproc))
if(nproc/2**N > 1) : N = N+1

eta01 = 1E20
eta02 = 5E20
eta03 = 1E21

Tm01 = 1690
Tm02 = 1710
Tm03 = 1750

k01 = 1E-11
k02 = 1E-9
k03 = 1E-7

rhocr1 = 3000
rhocr2 = 2800
rhocr3 = 2800

kcr1 = 3
kcr2 = 2.75
kcr3 = 2

A1 = 300E3
A2 = 200E3
A3 = 4

V1 = 3E-6
V2 = 6E-6
V3 = 5

fmag1 = 0.9
fmag2 = 0.5
fmag3 = 1

fbase1 = 0
fbase2 = 1
fbase3 = 0.5

CH2O1 = 0.02
CH2O2 = 0.06
CH2O3 = 0.05

dDl1 = 1E3
dDl2 = 3E3
dDl3 = 10E3


dTc01 = 100
dTc02 = 300
dTc03 = 200



table1 = [[eta01,Tm01,k01,rhocr1,kcr1,A1,V1,fmag1,fbase1,CH2O1,dDl1,dTc01]]
table2 = [eta02,Tm02,k02,rhocr2,kcr2,A2,V2,fmag2,fbase2,CH2O2,dDl2,dTc02]
table3 = [eta03,Tm03,k03,rhocr3,kcr3,A3,V3,fmag3,fbase3,CH2O3,dDl3,dTc03]
n_table = len(table2)

table1_bis = np.copy(table1)
table1_bis[0][0] = table2[0]
TABLE = np.append(table1,table1_bis,0)


for i in range(1,N) :
    TABLE_BIS = np.copy(TABLE)
    TABLE_BIS[:,i] = table2[i]
    TABLE = np.append(TABLE,TABLE_BIS,0)


# for i in range(2**N) :
#     print(i+1)
#     print(TABLE[i,N-1])
    

TABLE_write = TABLE[0:nproc,:]
#print(np.shape(TABLE_write))

np.savetxt("MPI_MCMC_These_entry.txt",TABLE_write,delimiter=' ', fmt = '%e')