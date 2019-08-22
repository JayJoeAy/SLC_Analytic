import numpy as np
def converge1(N,M,Vm,Vm_NEW,Vel_tol,V_iter,X):
    print("Converge1")
    V_iter = V_iter + 1
    V_flag  = 0
    
    for  i in range(1,M):
        for j in range(0,N):
            X[j,i]= np.abs(Vm[j,i]-Vm_NEW[j,i])
    
    for i in range(1,M):
        for j in range(0,N):
            if X[j,i] > Vel_tol:
                V_flag = 1
    
    if V_flag == 1:
        Vm = Vm_NEW
    
    return Vm,V_flag,V_iter
