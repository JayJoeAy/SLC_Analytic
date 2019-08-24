# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 19:40:26 2019

@author: Erfan
"""
#Profile Class
import numpy as np
class Mdot():
    def __init__(self,length,inRad,slope,r_h,M,N,rho_in,Vz_inT):
        mdot_init=self.SetMdot(rho_in,Vz_inT,N,inRad,r_h,length,M)
        self.length,self.inRad,self.slope,self.r_h,self.M,self.N=length,inRad,slope,r_h,M,N
        
    def SetMdot(self,rho_in,Vz_inT,N,inRad,r_h,length,M):
        self.mdot_tot=rho_in[0,0]*np.pi*(inRad**2-r_h**2)*Vz_inT
        self.mdot_per=self.mdot_tot/(N-1)
        self.z=np.linspace(0,length,M)
        self.r=np.zeros((N,M),dtype=float)
        self.rConstant=self.mdot_per[0]/(rho_in[0,0]*Vz_inT*np.pi)
        
class Nozzle(Mdot):     
    def SetRad(self,):
        length,inRad,slope,r_h,M,N=self.length,self.inRad,self.slope,self.r_h,self.M,self.N
        r_c=slope*self.z+inRad
        r=np.zeros((N,M),dtype=float)
        r[N-1,:]=r_c
        r[0,:]=r_h
        for i in range(1,N-1) :
            r[i,0]=np.sqrt(self.rConstant+r[i-1,0]**2)
    #    
        for i in range(1,M):
            for j in range(1,N-1):
                r[j,i]=( (r[j,0]-r[0,0])/(r[N-1,0]-r[0,0]) ) * (r[N-1,i]-r[0,i])+r[0,i]
                r[N-1,:]=r_c
            r[0,:]=r_h
        return r,self.mdot_tot,self.mdot_per