# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:49:28 2019

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt 
from scipy import optimize 
import myutil
out = myutil.read_fresco_res('02_BCCresults.dat')
for key in out.keys():
    out[key]=np.array(out[key])

free = {}
for n in [5,6,7,8,9,10,11]:
  free[n]= out[0][n-5,:]; #free[L][i] i=3 bcc_tr i=4 thermal 

out[0]; # free[i,j] i= L index, j = method index    
LL={}
for n in [5,6,7,8,9,10,11]:
  LL[n]=out[n-4];

#  B_tr
B_tr={}
for n in [5,6,7,8,9,10,11]:
  B_tr[n]=LL[n][:,3:]/free[n][3]

# B_fg
B_fg={}
for n in [5,6,7,8,9,10,11]:
  B_fg[n]=LL[n][:,3:]/free[n][4]

# 
Lt = np.arange(100,700,100)

def ffunc(t,E0,b,eta):
    return E0+b*np.exp(-eta*t)

p0list={}
p0list[5]=[0.35,0.1,0.01]
p0list[6]=[0.35,0.1,0.01]
p0list[7]=[0.35,0.1,0.01]
p0list[8]=[0.4,0.01,0.001]
p0list[9]=[0.35,0.01,0.001]
p0list[10]=[0.36,0.1,0.1]
p0list[11]=[0.36,0.1,0.1]
for n in [5,6,7,8,9]:
  #params, params_covariance = optimize.curve_fit(ffunc,Lt,B_tr[n][:,0],p0=p0list[n],sigma=B_tr[n][:,1])
  params, params_covariance = optimize.curve_fit(ffunc,Lt,B_tr[n][:,0],p0=p0list[n])
  print(params)
  print(np.sqrt(params_covariance[0,0]))
  xx=np.arange(0,700,10);
  ff=ffunc(xx,*params);
  plt.errorbar(Lt ,B_tr[n][:,0],B_tr[n][:,1],fmt='.',label='tr L=%i'%(n))
  plt.plot(xx ,ff)
plt.legend()
plt.title(r'$\beta^{finite}_{33,33}$')
plt.xlabel('Lt (a_t unit)')
plt.ylabel(r'$\beta^{finite}_{33,33}$')

for n in [5,6,7,8,9,10,11]:
  params, params_covariance = optimize.curve_fit(ffunc,Lt,B_fg[n][:,0],p0=p0list[n])
  print(params)
  xx=np.arange(0,700,10);
  ff=ffunc(xx,*params);
  plt.errorbar(Lt ,B_fg[n][:,0],B_fg[n][:,1],fmt='.',label='tr L=%i'%(n))
  plt.plot(xx ,ff)
plt.legend()
plt.title(r'$\beta^{thermal}_{33,33}$')
plt.xlabel('Lt (a_t unit)')
plt.ylabel(r'$\beta^{thermal}_{33,33}$')    
    
    
    





