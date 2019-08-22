
from scipy.interpolate import interp1d
def radial(i,N,M,mass_frac,r,Vm,gamma,R,T,delta_z,x,y,VM_NEW_i,r_new_i):
    
#    x=np.zeros((N))
#    y=np.zeros((N))
#    VM_NEW_i=np.zeros((N,M))
    
    
    mass_frac[N-1,0:M-1]=1
    mass_frac[0,0:M-1]=0
    for j in range(0,N):
        y[j]=r[j,i]
        x[j]=mass_frac[j,i]
        if x[j]>=1:
            x[j]=1
#    mass_frac_interp=mass_frac(1:N-2,1
    f1=interp1d(x,y)
    r_new_i[1:N-1,i] = f1(mass_frac[1:N-1,0])
    
    tsum=0
    
    for j in range(0,N):
        Mach2_i=(Vm[j,i]**2)  / (gamma * R * T[j,i])
        tsum = tsum + Mach2_i
    
    Ave_mach_i = tsum / N
    
#    del tsum
    RFr_i = 1/(1+ (5.0/24) * (1.0-Ave_mach_i**2) * ((r[N-1,i]-r[0,i]) / delta_z))
    
    for j in range(1,N-1):
        r_new_i[j,i]= r[j,i]+RFr_i*(r_new_i[j,i]-r[j,i])
    
    r_new_i[0,i]=r[0,i]
    r_new_i[N-1,i]=r[N-1,i]
    
    for j in range(0,N):
        y[j] = Vm[j,i]
        x[j] = r[j,i]
    
    f2=interp1d(x,y)
    VM_NEW_i[0:N-1,i]=f2(r_new_i[0:N-1,i])
    
    return r_new_i,VM_NEW_i,RFr_i,Ave_mach_i
    