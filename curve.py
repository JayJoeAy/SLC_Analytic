import numpy as np
from deriv import deriv
def curve(i,r,z,N,M,C_i,phi_i,dr2dz2,drdz,gamma,delta_z):
    if i==M-1:
        for j in range(0,N):
            drdz[j,M-1] = deriv(i,r[j,0:M],z,M,1)
            dr2dz2[j,M-1] = deriv(i,r[j,0:M],z,M,2)
            C_i[j,M-1] = - dr2dz2[j,M-1] / ((1 + drdz[j,M-1] ** 2) ** 1.5)
    elif i==0:
        for j in range(0,N):
            drdz[j,0] = deriv(i,r[j,0:M],z,M,1)
            dr2dz2[j,0] = deriv(i,r[j,0:M],z,M,2)
            C_i[j,0] = - dr2dz2[j,0] / ((1 + drdz[j,0] ** 2) ** 1.5)
    else :
        for j in range (0,N):
            dr2dz2[j,i] = deriv(i,r[j,0:M],z,M,2)
            drdz[j,i] = deriv(i,r[j,0:M],z,M,1)
            C_i[j,i] = - dr2dz2[j,i] / ((1 + drdz[j,i] ** 2) ** 1.5)            
        
    for j in range(0,N):
        phi_i[j,i]=np.arctan(drdz[j,i])
        if j==0:
            gamma[j,i]=np.arctan(delta_z/(r[j+1,i]-r[j,i]))
        else:
            gamma[j,i]=np.arctan(delta_z/(r[j,i]-r[j-1,i]))
        
    return C_i,phi_i,gamma
