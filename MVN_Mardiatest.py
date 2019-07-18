### Mardia normality test
import numpy as np
from scipy as stats

def Mardiatest(n, k, Data):
    def gamma(n, Data):
        S = np.cov(Data.T, bias = True)
        gamma_1 = 0; gamma_2 = 0; m_ij = 0
        for i in range(n):
            for j in range(n):
                m_ij = np.matmul(np.matmul((Data[i]- np.mean(Data, axis = 0)).transpose(), inv(S)), (Data[j]-np.mean(Data, axis = 0)))
                gamma_1 += m_ij
                if j == i:
                    gamma_2 += pow(m_ij, 2)

        #print(beta, A, B, C, T)
        return gamma_1/pow(n,2), gamma_2/n      
    gamma1, gamma2 = gamma(n, Data)
    gamma1_df = k*(k+1)*(k+2)/6
    gamma2_mean = k*(k+2)
    gamma2_var = 8*k*(k+2)/6
    pvalue_gamma1 = 1- stats.chi2.cdf(gamma1*n/6, gamma1_df)
    pvalue_gamma2 = 1- stats.norm.cdf(gamma2, loc = gamma2_mean, scale = gamma2_var)
    return np.minimum(pvalue_gamma1, pvalue_gamma2)
    