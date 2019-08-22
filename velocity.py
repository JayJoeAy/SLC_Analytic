import numpy as np
def velocity(i,SL_mid,Vm,phi,C,N,dvmdm,dHodq,dsdq,
    drvtdr,r,RFv,rho,mdot_tot,mass_flow_tol,mass_iter,
    FM,FN,T_0,M,Old_Vm2,VM2NEW,psi,A,B,Tube_Flow,mass_frac,A_Bar):
    
#    M=80
#    Old_Vm2=np.zeros((N))
#    VM2NEW=np.zeros((N))
#    epsilon=np.zeros((N,M))
#    psi=np.zeros((N,M))
#    A=np.zeros((N))
#    B=np.zeros((N))
#    Tube_Flow=np.zeros((N))
#    FN=np.zeros((N,M))
#    FM=np.zeros((N,M))
#    mass_frac=np.zeros((N,M))
#    A_Bar=np.zeros((N)
    flag_mass=1
    Vm_MID = Vm[SL_mid,i]
    
    for j in range(0,N):
        Old_Vm2[j]=Vm[j,i]**2
    
    while flag_mass == 1:
        VM2NEW[SL_mid]=Vm_MID**2
        
        for j in range (0,N):
            A[j]= 2 * (np.cos(phi[j,i] + psi[j,i]) * C[j,i] - np.sin(phi[j,i] +
            psi[j,i]) / Vm[j,i] * dvmdm[j,i])
            B[j] = 2 * (dHodq[j,i] - T_0[j,i] * dsdq[j,i] - 1 /
            (2 * r[j,i] ** 2) * drvtdr[j,i] - (np.sin(phi[j,i] + psi[j,i]) * 
            FM[j,i] + np.cos(phi[j,i] + psi[j,i]) * FN[j,i])/rho[j,i])
        
        for j in range(SL_mid+1,N):
            A_Bar[j] = (A[j-1] + A[j]) / 2
            VM2NEW[j] = np.exp( -A_Bar[j] * (r[j,i] - r[j-1,i]) ) * (VM2NEW[j-1] + ((r[j,i] - r[j-1,i]) * B[j-1] / 2)) + B[j] * (r[j,i] - r[j-1,i]) / 2
            
        for j in range(SL_mid-1,-1,-1):
            A_Bar[j] = (A[j+1] + A[j]) / 2
            VM2NEW[j] = np.exp( A_Bar[j] * (r[j+1,i] - r[j,i]) ) * (VM2NEW[j+1] - ((r[j+1,i] - r[j,i]) * B[j+1] / 2)) - B[j] * (r[j+1,i] - r[j,i]) / 2
        
        for j in range(0,N):
            VM2NEW[j] = Old_Vm2[j] + RFv * (VM2NEW[j]-Old_Vm2[j])
            Old_Vm2[j] = VM2NEW[j]
            
        for j in range(0,N):
            if VM2NEW[j]<0:
                VM2NEW[j]=1.0
            
            Vm[j,i]=np.sqrt(VM2NEW[j])
        
        Mass_Flow=0
        Tube_Flow[1]=0
        
        for j in range(1,N):
            Tube_Flow[j] = (rho[j-1,i] + rho[j,i]) / 2 * np.pi * ((r[j,i]) ** 2 -
            (r[j-1,i]) ** 2) * (Vm[j,i] * np.cos(phi[j,i]) + Vm[j-1,i] *
            np.cos(phi[j-1,i])) / 2
            Mass_Flow = Mass_Flow + Tube_Flow[j]
        
        if np.abs(Mass_Flow-mdot_tot[0]) < mass_flow_tol :
            flag_mass =0
            mass_frac[0,i]=0
            mass_frac[N-1,i]=1
            
            for j in range(1,N-1):
                mass_frac[j,i]=Tube_Flow[j]/Mass_Flow+ mass_frac[j-1,i]
        if Mass_Flow/mdot_tot[0] < 1:
            Vm_MID= np.minimum( 1.2 * Vm_MID,(2-(Mass_Flow/mdot_tot[0]))*Vm_MID )     
        if Mass_Flow/ mdot_tot[0] > 1 :
            Vm_MID= np.maximum( 0.8 * Vm_MID,(2-(Mass_Flow/mdot_tot[0]))*Vm_MID )
        mass_iter=mass_iter+1
    return Vm,mass_frac,mass_iter,Mass_Flow,A,B
