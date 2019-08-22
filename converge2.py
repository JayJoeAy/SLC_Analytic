import numpy as np
def converge2(r,r_new,r_tol,M,N,VM_NEW,Vm,r_iter):
    R_flag = 0 
    r_iter = r_iter +1
    X=np.zeros((N,M))
    
    for i in range (1,M):
        for j in range (0,N):
            X[j,i]= np.abs(r[j,i]-r_new[j,i])
    
    for i in range (1,M):
        for j in range (0,N):
            if X[j,i] > r_tol:
                R_flag = 1 
    
    if R_flag == 1 :
        for i in range (1,M):
            for j in range (0,N):
                Vm[j,i] = VM_NEW[j,i]
                r[j,i] = r_new[j,i]
        
    else:
        Vm=Vm
        r=r
    
    return R_flag,Vm,r,r_iter
    
    
    