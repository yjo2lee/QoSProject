# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

### multivariate normal test using Henze-Zirkler method
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.stats import norm

def Tvariable(n, k, Data):
    S = np.cov(Data.T, bias = True)
    #beta = 0.5
    beta = (1/np.sqrt(2))*((n*(2*k+1)/4)**(1/(k+4)))
    C = n*(1 + 2*beta**2)**(-k/2)  
    B = 0; T = 0
    for i in range(n):
        A = 0
        for j in range(n):
            A += np.exp(-(beta**2/2)*(np.matmul(np.matmul((Data[i]- Data[j]).transpose(), inv(S)), (Data[i]-Data[j]))))
        B = np.exp(-(beta**2/(2*(1+beta**2)))*(np.matmul(np.matmul((Data[i]- np.mean(Data, axis = 0)).transpose(), inv(S)), (Data[i]-np.mean(Data, axis = 0)))))
        T += (1/n)*A -2*(1 + beta**2)**(-k/2) * B
    #print(beta, A, B, C, T)
    return beta, T + C
        

#def Moments(beta, k):
#   w = (1+pow(beta, 2))*(1+3*pow(beta, 2))
#   a = 1+2*pow(beta, 2)
#   Expec_T = 1- pow(a,-k/2)*(1 + k*pow(beta, 2/a) + k*(k+2)*pow(beta,4))/2/pow(a,2)
#   Var_T = 2*pow(1 + 4*pow(beta, 2),-k/2) + 2*pow(a,-k)*(1+ 2*k*pow(beta,4))/pow(a,2) + 3*k*(k+2)*pow(beta,8)/4/pow(a,4)-4*pow(w,(-k/2))*(1 + (3*k*pow(beta,4)/2/w) + (k*(k+2)*pow(beta,8)/2/pow(w,2)))
#   #print(w, Expec_T, Var_T)
#   return Expec_T, Var_T

def Moments(beta, k):
    w = (1+beta**2)*(1+3*beta**2)
    Expec_T = 1-((1 + 2*beta**2)**(-k/2)*(1 + k*beta**2/(1+2*beta**2) + k*(k+2)*beta**4/2/(1+2*beta**2)**2))
    Var_T = 2*(1 + 4*beta**2)**(-k/2) + 2*((1+2*beta**2)**(-k))*(1+ 2*k*beta**4/(1+2*beta**2)**2 + (3*k*(k+2)*beta**8)/(4*(1+2*beta**2)**4))-(4*w**(-k/2))*(1 + (3*k*beta**4/2/w) + (k*(k+2)*beta**8/2/w**2))
    return Expec_T, Var_T


def jointpval(n, k, Data):
    beta, T = Tvariable(n, k, Data)
    Expec_T, Var_T = Moments(beta, k)
    VZ = np.log(1+ Var_T/(Expec_T**2))
    SD = np.sqrt(VZ)
    EZ = np.log(np.sqrt(Expec_T**4/(Var_T+Expec_T**2)))
    stand_z = (T - EZ)/(np.sqrt(VZ))
    print(T, Expec_T, Var_T, EZ, SD, stand_z)    
    return joint_pval


Y = np.random.normal(0,1, (50,2))
X = np.array([[-7,5],[-3,15]])
A = np.matmul(X,Y.T).T
pval_test = jointpval(50, 2, A)

print(pval_test)
