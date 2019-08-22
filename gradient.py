import numpy as np
from deriv import deriv
def gradient(i,r,T_0,N,Vm,rho,V_t,s,H_0,Cp,R,gamma,T,SL_mid,p_0,V_tot,P,delta_m,rV_t2,M,omega_m,ct,ch,Old_Vm,dsdq,dHodq,drvtdr,dvmdm,delta_po,rel_f,delta_s,H):

    "ct and ch are taken to be zero and three first for's are neglected"
    
#    for j in range(0,N):
#        if i==0:
#            delta_s[j,i]=Cp*((gamma/(gamma-1))*
#            np.log((p_0[j,i]-np.abs(delta_po[j,i]))/p_0[j,i])-
#            np.log(T_0[j,i]/T_0[j,i]))
#        else:
#            delta_s[j,i]=Cp*((gamma/(gamma-1))*
#            np.log((p_0[:,i-1]-np.abs(delta_po[j,i]))/p_0[j,i-1])-
#            np.log(T_0[j,i]/T_0[j,i-1]))
    delta_s[:,0]=Cp*((gamma/(gamma-1))*np.log((p_0[:,0]-np.abs(delta_po[:,0]))/p_0[:,0])-np.log(T_0[:,0]/T_0[:,0]))        
    delta_s[:,1:]=Cp*((gamma/(gamma-1))*np.log((p_0[:,:-1]-np.abs(delta_po[:,1:]))/p_0[:,:-1])-np.log(T_0[:,1:]/T_0[:,:-1]))
        
#    for j in range(0,N):
#        if i==0:
#            s[j,i]=s[j,i]-delta_s[j,i]
#        else:
#            s[j,i]=s[j,i-1]-delta_s[j,i]
    s[:,0]=s[:,0]-delta_s[:,0]
    s[:,1:]=s[:,:-1]-delta_s[:,1:]
            
    for j in range (0,N):
        dsdq[j,i]=deriv(j,s[0:N,i],r[0:N,i], N ,1)
        
    "Updating all other flow field variables according to the new pressure and entropies"
    V_tot=np.sqrt(Vm**2+V_t**2)
    H=H_0-V_tot**2/2
    T=H/Cp
    P=(T**(Cp/R))*np.exp(-s/R)
    p_0=P*(T_0/T)**(gamma/(gamma-1))
    rho=(P)/(R*T)
        
    for j in range (0,N):
        dHodq[j,i]=deriv(j,H_0[0:N,i],r[0:N,i],N,1)
        
    for j in range (0,N):
        dsdq[j,i]=deriv(j,s[0:N,i],r[0:N,i],N,1)
        
    for j in range (0,N):
        drvtdr[j,i]=deriv(j,rV_t2[0:N,i],r[0:N,i],N,1) 
    
    dvmdm[:,1:-1] = (delta_m[:,2:] / (delta_m[:,1:-1] + delta_m[:,2:])) * (Vm[:,2:] - Vm[:,1:-1]) /delta_m[:,2:] + (delta_m[:,2:] / (delta_m[:,1:-1] + delta_m[:,2:])) * (Vm[:,1:-1] - Vm[:,:-2]) / delta_m[:,1:-1]
    dvmdm[:,0]=(delta_m[:,0] + delta_m[:,1]) / delta_m[:,0] *(Vm[:,1]-Vm[:,0])/delta_m[:,1] -(delta_m[:,0] / (delta_m[:,1] + delta_m[:,0])) * (Vm[:,2] - Vm[:,0]) / delta_m[:,1]
    dvmdm[:,-1] = (delta_m[:,-1] + delta_m[:,-2]) / delta_m[:,-1] *(Vm[:,-1]-Vm[:,-2])/delta_m[:,-1] -(delta_m[:,-1] / (delta_m[:,-1] + delta_m[:,-2])) * (Vm[:,-1] -Vm[:,-3]) / delta_m[:,-2]
    
    return rV_t2,s,V_t,T_0,H,T,P,p_0,V_tot,rho,dsdq,drvtdr,dHodq,dvmdm,delta_s,delta_po
