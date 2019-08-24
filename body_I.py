from PyQt5.uic import loadUiType
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QMainWindow, QTextEdit, QAction, QFileDialog, QTableWidget,QTableWidgetItem

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure

from bokeh.plotting import figure, output_file, show
from bokeh.io import output_file

output_file('outfile.htm')

import numpy as np
import matplotlib.pyplot as plt  
from deriv import deriv
from curve import curve
from converge1 import converge1
from gradient import gradient
from velocity import velocity
from radial import radial
from converge2 import converge2

Ui_MainWindow, QMainWindow = loadUiType('window.ui')

class Main(Ui_MainWindow,QMainWindow):
    def __init__(self, ):
        super(Main,self).__init__()
        self.setupUi(self)
        fig=Figure()
        self.addmpl(fig)
        
    def addmpl(self,fig):
        self.canvas=FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar=NavigationToolbar(self.canvas,self,coordinates=True)
        self.mplvl.addWidget(self.toolbar)
    def rmmpl(self,):
        fig1=Figure()
        ax1=fig1.add_subplot(111)
#        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
#        self.mplvl.removeWidget(self.toolbar)        
        self.toolbar.close()

def SlcMain():

    "Velocity relaxation factor"
    RFv=0.15
    mass_flow_tol = 1e-2
    Vel_tol=0.0001
    r_tol=1e-2
    
    z=np.linspace(0,2,50)
    M=len(z)
    N=50
    #T_0_in=288.7 #Kelvin
    T_0_in=float(main.l_t.text())
    
    #P_0_in=101352.932 #Pascal
#    P_0_in=2116.8
    P_0_in=float(main.l_p.text())
    
    Vr_in=0.0
    Vz_inT=float(main.l_v.text())
    gamma=1.4
    omega_m=0.0
    omega1=2
    ct=0
    ch=0
    
    "Casing and hub profile"
    r_m=1
    r_h=np.zeros((M))+0.4
    a=0.2
    b=2.5
    d=1
#    r_c=r_m+0.25*(d-z)
    r_c=r_m-a*np.tanh(b*(z-d))
    #r_c=np.zeros((M))+1
    #r_c=a*z+r_h+0.2
    #plt.plot(z,r_c)
    
    "Angles"
    epsilon=np.zeros((N,M))
    psi=np.zeros((N,M))
    
    Cp=6009.07
    
    V_in = np.sqrt(Vz_inT**2+Vr_in**2)
    T_in=np.zeros((N,M),dtype=float)
    T_in[0:,0]=T_0_in-(V_in**2)/(2*Cp)
    
    R=1715.84
    #R=cp.PropsSI('gas_constant','T',T_0_in,'P',P_0_in,AIR)
    
    #s_in1=cp.PropsSI('S','T',T_0_in,'P',P_0_in,AIR)
    #s_ref=cp.PropsSI('S','T',1,'P',1,AIR)
    s_in=Cp*np.log(T_0_in)-R * np.log(P_0_in)
    
    
    #rho_in=cp.PropsSI('D','T',T_0_in,'P',P_0_in,AIR)
    rho_in=(T_in**((Cp-R)/(R)) * np.exp(-s_in/R))/R
    
    mdot_tot=rho_in[0,0]*np.pi*(r_c[0]**2-r_h**2)*Vz_inT
    mdot_per=mdot_tot/(N-1)
    
    #****************************حدس های اولیه******************************************
    
    ###############################################################################
    
    #حدس اولیه ای برای شعاع با توجه هب پرفایل کیس و هاب
    
    ###############################################################################
    r=np.zeros((N,M),dtype=float)
    r[N-1,:]=r_c
    r[0,:]=r_h
    for i in range(1,N-1) :
        r[i,0]=np.sqrt(mdot_per[0]/(rho_in[0,0]*Vz_inT*np.pi)+r[i-1,0]**2)
    
    for i in range(1,M):
        for j in range(1,N-1):
            r[j,i]=( (r[j,0]-r[0,0])/(r[N-1,0]-r[0,0]) ) * (r[N-1,i]-r[0,i])+r[0,i]
    
    delta_z=z[1]-z[0]
    
    SL_mid=int(N/2)
    
    ###############################################################################
    
    # مقادیر سرعت V_theta ورودی 
    # همان پارامتر whirl»
    
    ###############################################################################
    #V_t=np.zeros((N,M))+50*r
    V_t=np.zeros((N,M))
    V_t[:,(0,1)]=50*r[:,(0,1)]
    V_t[:,(4,5)]=50*r[:,(4,5)]
    #whirl=(r*V_t)**2
    #-----------------------------------------
    
    Vz_in=np.zeros((N,M))
    Vz_in[:,0]=Vz_inT
    Vz=np.zeros((N,M))
    Vz=Vz+Vz_in
    U=np.zeros((N,2*M))
    #U[:,81:]=50*r
    #"Plotting the streamlines"
    #plt.subplot(1, 2, 1)
    #for i in range(0,N):
    #    plt.plot(z,r[i,:])
    #plt.savefig('first one')
    Vz_in=Vz_in[:,0]
    
    ######################################################################################
    #  این قسمت از معادلات تنها برای وقتی است که ما برای سرعت ورودی مقدار V_theta در#
    # نظر بگیریم در این صورت است که مقدار V_z باید تغییر کند وگرنه که همان مقدار ثابت باید #
    # باقی بماند. هدف استفاده هم این است که از معادله تعادل شعاعی ساده شده برای     #
     # بدست آوردن سرعت محوری با توجه به مقدار سرعت V_t ورودی برای محاسبه مقادیر    #
     # آنتالپی و دمای سکون فقط ورودی است.خط 91                                        #
    # در واقع پارامتر rv_t^2 که همان whirl باشد از تابع مشتق گیری بدست می آید          #
    ######################################################################################
    "Simplified radial eq"
    "Below midline"
    for j in range(SL_mid-1,-1,-1):
        der1=deriv(j,(r[:,0]** 2 * V_t[:,0]** 2),r[:,0],N,1)
        Vz_in[j] =np.sqrt(Vz_in[j+1] ** 2 - 0.5 * (r[j,0] - r[j+1,0]) *
        (1 / (r[j,0] ** 2) * der1
        + 1 / (r[j+1,0] ** 2) * deriv(j+1,(r[:,0] ** 2 *
        V_t[:,0] ** 2),r[:,0],N,1)))
    
    del j
    "Above midline"
    for j in range(SL_mid+1,N):
        der1=deriv(j,(r[:,0] ** 2 * V_t[:,0] ** 2),r[:,0],N,1)
        Vz_in[j] = np.sqrt(Vz_in[j-1] ** 2 - 0.5 * (r[j,0] - r[j-1,0]) * 
        (1 / (r[j,0] ** 2) * der1
        + 1 / (r[j-1,0] ** 2) * deriv(j-1,(r[:,0] ** 2 *
        V_t[:,0] ** 2),r[:,0],N,1)));
    #-------------------------------------------------------------------------------------
    
    ###############################################################################
    
    #با توجه به قضیه استفاده از نسبت های دبی یک حدس اولیه از سرعت محوری بدست می آید
        
    ###############################################################################
    "Find the proportion of the mass flow rate in each stream tube"
    mass_frac_temp=np.zeros((N,M))
    mass_frac = np.zeros((N,M))
    for j in range(1,N):
        mass_frac_temp[j]=mdot_per[1]/mdot_tot[1]+mass_frac_temp[j-1]    
    mass_frac[:,0]=mass_frac_temp[:,0]
    
    for i in range(1,M):
        Vz[0:N,i]=mdot_tot[0]/(rho_in[0,0]*np.pi*(r[N-1,i]**2-r[0,i]**2)) 
    #------------------------------------------------------------------------------
    
    ###############################################################################
    # این جا حدس اولیه برای تمام شبه نرمال های گرید بندی انجام میشود با این   
    # فرض که فشار و انتالپی و دمای سکون در باقی شبه نرمال ها برابر با فشار و 
    # انتالپی سکون ورودی است.                                             
    ###############################################################################
    " فشار سکون "
    p_0=np.zeros((N,M),dtype=float)+P_0_in
    H_0_in=np.zeros((N))+Cp*T_0_in
    H_in=H_0_in-0.5*Vz_in**2
    "دمای سکون"
    T_0=np.zeros((N,M))+T_0_in
    V_tot=np.sqrt(Vz**2+V_t**2)
    T=T_0-(V_tot**2)/(2*Cp)
    "انتالپی سکون"
    H_0=np.zeros((N,M))+H_0_in[1]
    H=Cp*T
    s=Cp*np.log(T_0)-R*np.log(p_0)
    P=np.exp((Cp/R*np.log(T))-s/R)
    rho=P/(R*T)
    #------------------------------------------------------------------------------
    
    #*********************************پایان حدس های اولیه*********************************
    
    
    mass_f=np.zeros((N,M))+mass_frac[0:N]
    #del T_0_in,T_in,H_0_in,H_in , P_0_in ,s_in, r_h, r_c, V_in, Vr_in, 
    mass_frac=mass_f
    #del mass_f
    
    "Setting mass iteration counter equal to zero"
    mass_iter=0
    "Setting velocity iteration counter equal to zero"
    V_iter=0
    "Setting Radial change iteration counter equal to zero"
    r_iter=0
    
    "Relative Components"
    W_t=np.zeros((N,2*M))
    
    "Setting the flags"
    V_flag=1
    R_flag=1
    
    C=np.zeros((N,M))
    C_i=np.zeros((N,M))
    phi_i=np.zeros((N,M))
    gamma_angle=np.zeros((N,M))
    dr2dz2=np.zeros((N,M))
    drdz=np.zeros((N,M))
    
    
    Old_Vm=np.zeros((N,M))
    delta_m=np.zeros((N,M))
    rV_t2=np.zeros((N,M))
    dsdq=np.zeros((N,M)) 
    dHodq=np.zeros((N,M)) 
    drvtdr = np.zeros((N,M)) 
    dvmdm = np.zeros((N,M)) 
    delta_po = np.zeros((N,M)) 
    rel_f =np.zeros((N,M))
    num=np.zeros((N,M))
    d1dm=np.zeros((N,M))
    delta_s=np.zeros((N,M))
    
    FN=np.zeros((N,M))
    FM=np.zeros((N,M))
    
    Old_Vm2=np.zeros((N))
    VM2NEW=np.zeros((N))
    A_Bar=np.zeros((N))
    A=np.zeros((N))
    B=np.zeros((N))
    Tube_Flow=np.zeros((N))
    
    X=np.zeros((N,M))
    
    x=np.zeros((N))
    y=np.zeros((N))
    VM_NEW_i=np.zeros((N,M))
    r_new_i=np.zeros((N,M))
        
    while R_flag==1:
        while V_flag==1:
            for i in range(0,M):
                "Calcluation of curvature and flow angles"
                C,phi,gamma_angle=curve(i,r,z,N,M,C_i,phi_i,dr2dz2,drdz,gamma_angle,delta_z)
            Vm=Vz/np.cos(phi)
            Vr=Vm*np.sin(phi)
            if mass_iter > 0:
    #            Vm[0:N,0]=Old_Vm[0:N,0]
                Old_Vm=Vm.copy()
            "Gradient"
    #        for i in range (0,M):
    #            for j in range (0,N):
    #                if i==M-1:
    #                    delta_m[j,M-1]=delta_m[j,M-2]
    #                    rV_t2[j,i]=(r[j,i]*V_t[j,i])**2
    #                else:
    #                    delta_m[j,i]=np.sqrt(((r[j,i+1]-r[j,i]))**2+
    #                    (z[i+1]-z[i])**2)
    #                    rV_t2[j,i]=(r[j,i]*V_t[j,i])**2
    #        delta_m2=np.zeros_like(delta_m)
            delta_m[:,:-1]=np.sqrt( (r[:,1:]-r[:,:-1])**2 + (z[1:]-z[:-1])**2 )
            delta_m[:,-1]=delta_m[:,-2]
            rV_t2=(r*V_t)**2
            
    #        for i in range (0,M):
            rV_t2,s,V_t,T_0,H,T,p,p_0,V_tot,rho,dsdq,drvtdr,dHodq,dvmdm,delta_s,delta_po=gradient(i,r,T_0,N,Vm,rho,V_t,s,H_0,Cp,R,gamma,T,SL_mid,p_0,V_tot,
            P,delta_m,rV_t2,M,omega_m,ct,ch,Old_Vm,dsdq,dHodq,drvtdr,dvmdm,delta_po,rel_f,delta_s,H)
            
            
            Old_Vm=Vm
    #        for i in range (1,M):
    #            mass_frac[0:N,i] = mass_frac[0:N,0]
    #            
    #        
            
            W_t=V_t
            
            TanBeta=W_t/Vm        
    #        for i in range(1,M):
    #            for j in range(0,N):
            num[:,1:]=r[:,1:]*Vm[:,1:]*TanBeta[:,1:]-r[:,:-1]*Vm[:,:-1]*TanBeta[:,:-1]
    #                delta_m[j,i]=np.sqrt( ((r[j,i]-r[j,i-1]))**2 +(z[i]-z[i-1])**2)
            d1dm=num/delta_m
            
    # In the Unbladed region blade forces are zero"
    #        for i in range(0,M):
    #            for j in range(0,N):
    #        Fthet=0
    #        FN=0
    #        FM=0
    #        
    #        for i in range(0,M):
    #            for j in range(0,N):
    #                Fthet=(Vm[j,i]/r[j,i])*d1dm
    #                FN[j,i]=Vm[j,i]/r[j,i] * (np.tan(rel_f[j,i])*np.tan(phi[j,i]+psi[j,i])+np.tan(epsilon[j,i])/np.cos(phi[j,i]+psi[j,i]))*d1dm[j,i]
    #                FM[j,i]=Vm[j,i]*np.tan(rel_f[j,i])/r[j,i]*d1dm[j,i]
    #                    
            for i in range(1,M):
                Vm_NEW,mass_frac,mass_iter,Mass_Flow,A,B=velocity(i,SL_mid,Vm,phi,C,N,dvmdm,dHodq,dsdq,drvtdr,r,RFv,rho,
                mdot_tot,mass_flow_tol,mass_iter,FM,FN,T_0,M,Old_Vm2,VM2NEW,psi,A,B,Tube_Flow,mass_frac,A_Bar)
            Vm=Old_Vm
            Vm,V_flag,V_iter=converge1(N,M,Vm,Vm_NEW,Vel_tol,V_iter,X)
            Vz = Vm*np.cos(phi)
        
        for i in range(1,M):
            r_NEW,Vm_NEW,RFr,Ave_mach=radial(i,N,M,mass_frac,r,Vm,gamma,R,T,delta_z,x,y,VM_NEW_i,r_new_i)
            
        r_NEW[0:N,0]=r[0:N,0]
        
    
        
        R_flag,Vm,r,r_iter=converge2(r,r_NEW,r_tol,M,N,Vm_NEW,Vm,r_iter)
        Vz= Vm * np.cos(phi)
        V_tot = np.sqrt(Vm**2 + V_t**2)
        T = T_0 - (V_tot **2)/(2*Cp)
        H = Cp * T_0 - ( V_tot **  2 ) / (2 * Cp)
        p = np.exp( (Cp / R * np.log(T)) - s / R)
        rho = p/ (R* T)
    
    Vz= Vm * np.cos(phi)
    Vr= Vm * np.sin(phi)
    V_tot = np.sqrt(Vm**2 + V_t**2)
    T = T_0 - (V_tot **2)/(2*Cp)
    H = Cp * T_0 - ( V_tot **  2 ) / (2 * Cp)
    Hst=Cp*T_0
    p = np.exp( (Cp / R * np.log(T)) - s / R)
    rho = p/ (R* T)
    
    #****************************************************** Bladed region ********************************** 
    "Casing and hub profile"
    r_m=1
    r_h=np.zeros((M))+0.4
    a=0.2
    b=2.5
    d=1
    r_c=r_m+0.25*(d-z)
#    TP=[("Vm", "@Vm")]
    p=figure()
    fig1=Figure()
    ax1=fig1.add_subplot(111)
    #plt.subplot(1, 2, 2)
    for i in range(0,N):
        p.line(z,r_NEW[i,:])
        ax1.plot(z,r_NEW[i,:])
    z_plot=np.ones((N,M))*z
    for u in range(0,M):
        ax1.plot(z_plot[:,u],r_NEW[:,u])
        p.line(z_plot[:,u],r_NEW[:,u])
    main.rmmpl()
    main.addmpl(fig1)
    
#    hover.tooltips=["Vm", "@Vm"]
#    show(p)
if __name__=='__main__':
    import sys
    from PyQt5 import QtGui
    
    app=QApplication(sys.argv)
    main=Main()
    
    main.show()
    
    main.actionNozzle.triggered.connect(lambda : NozzleProfile())
    main.actionDiffuser.triggered.connect(lambda : DiffuserProfile())
    main.actionCNV_DIV.triggered.connect(lambda : CNV_DIVProfile())
    
    main.ButStart.clicked.connect(lambda : SlcMain())
    main.ExitBut.clicked.connect(QApplication.instance().quit)
    sys.exit(app.exec_())
        