# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 20:18:14 2023

@author: 86131
"""

import numpy as np
import matplotlib.pyplot as plt


#1.生成矩阵D与F，并解方程，得到系数c
def DF(N):
    D=np.zeros((N,N))
    F=np.zeros((N,))
    for i in range(1,len(F)+1):
        F[i-1]=1/((3+i)*(4+i))
        for j in range(i,len(F)+1):
            D[i-1][j-1]=D[j-1][i-1]=(2*i*j)/((i+j+1)*(i+j)*(i+j-1))
            
    c=np.linalg.inv(D)@F
    return c

# #生成u的近似函数u_N
# def uN(c,N,x):
#     uN=0
#     for i in range(1,N+1):
#         uN=uN+c[i-1]*(x**i * (1-x))
#     return uN

#2.计算f
def f(N,c,x):
    uN=0
    for i in range(1,N+1):
        uN=uN+c[i-1]*(x**i * (1-x))
    f=((x*(1-x**3))/12 - uN)**2
    return f

#3.计算误差(估计两函数差的平方之积分)
#k为最大二分次数
def e(N,c,k):
    h=1/2**k
    T_0_0=(h/2)*(f(N,c,0)+f(N,c,1))
    for i in range(1,2**k):
        T_0_0=T_0_0+h*f(N,c,i/2**k)
    T_0_1=0.5*T_0_0
    for i in range(2**k):
        T_0_1=T_0_1+(h/2)*f(N,c,h*i+h/2)
    T_1_k=(4/3)*T_0_1-(1/3)*T_0_0
    e=T_1_k*108
    return e

           
#4.将N不断增大，观察走势
N=10
err=np.zeros((N,))
for n in range(1,N+1):
    c_0=DF(n)
    err[n-1]=e(n,c_0,5)

plt.figure()
plt.plot(np.arange(1,N+1),err)

plt.figure()
plt.plot(np.arange(3,N+1),err[2:])

c=DF(3)
















    


