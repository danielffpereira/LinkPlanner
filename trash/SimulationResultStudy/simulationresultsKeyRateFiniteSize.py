import numpy as np
import matplotlib.pyplot as plt

Xi=[0.1,0.11]


def g(x):
    return (x+1)*np.log2(x+1)-x*np.log2(x)

d=np.arange(0,30.01,.001)
eta=1
T=eta*10**(-.02*d)
Va=1;

DS11=np.array([0,10.1,15.1,5.1])
KS11=np.array([.0148,.00393,.000462,.00815,])
KS11min=np.array([.0152,-.0474,-.0713,-.0200])
KS11max=np.array([-.00514,.0424,.0645,.0199])

DS10=np.array([0,10,20,25,26,5,15])
KS10=np.array([.0283,.0143,.00496,.00167,.000293,.0195,.00814])
KS10min=np.array([.0289,-.0359,-.0844,-.104])
KS10max=np.array([0.00801,.0515,.0951,.129])

zepe=6.4670
SIZE=10**10

z=.5*np.exp(-Va/2)*np.array([np.cosh(Va/2)+np.cos(Va/2),np.sinh(Va/2)+np.sin(Va/2),np.cosh(Va/2)-np.cos(Va/2),np.sinh(Va/2)-np.sin(Va/2)])
Z=Va*(z[3]**(3/2)/z[0]**.5+z[0]**(3/2)/z[1]**.5+z[1]**(3/2)/z[2]**.5+z[2]**(3/2)/z[3]**.5)
sigmaZ=np.array([[1,0],[0,-1]])
S=np.zeros([len(T)])
Smin=np.zeros([len(T)])
for j in range(0,len(Xi)):
    xi=Xi[j];
    tmin=np.sqrt(T/2)-zepe*np.sqrt((1+T/2*xi)/(SIZE*Va))
    sigmamax=1+T/2*xi+zepe*(1+T/2*xi)*np.sqrt(2)/np.sqrt(SIZE)
    
    
    SNR=T/2*Va/(1+T/2*xi)
    Iab=np.log2(1+SNR)
    
    SNRmin=tmin**2*Va/sigmamax
    Iabmin=np.log2(1+SNRmin)
    
    for i in range(0,len(T)):
        gammaA=(1+Va)*np.eye(2)
        gammaB=(T[i]/2*Va+1+T[i]/2*xi)*np.eye(2)
        gammaC=np.sqrt(T[i]/2)*Z*sigmaZ
        gamma=np.vstack([np.hstack([gammaA, gammaC]), np.hstack([gammaC, gammaB])])
        
        gammaAB_B=(Va+1-T[i]/2*Z**2/(T[i]/2*Va+2+T[i]/2*xi))*np.eye(2)
        
        
        gammaAmin=(1+Va)*np.eye(2)
        gammaBmin=(tmin[i]**2*Va+sigmamax[i])*np.eye(2)
        gammaCmin=tmin[i]*Z*sigmaZ
        
        gammamin=np.vstack([np.hstack([gammaAmin, gammaCmin]), np.hstack([gammaCmin, gammaBmin])])
        gammaAB_Bmin=(Va+1-tmin[i]**2*Z**2/(tmin[i]*Va+1+sigmamax[i]))*np.eye(2)
        
        
        Delta=np.linalg.det(gammaA)+np.linalg.det(gammaB)+2*np.linalg.det(gammaC)
        mu1=(.5*(Delta+(Delta**2-4*np.linalg.det(gamma))**.5))**.5
        mu2=(.5*(Delta-(Delta**2-4*np.linalg.det(gamma))**.5))**.5
        k=np.linalg.eig(gammaAB_B)
        mu3=k[0]
        mu3=mu3[0]
        n1=(mu1-1)/2
        n2=(mu2-1)/2
        n3=(mu3-1)/2
        S[i]=g(n1)+g(n2)-g(n3)
        
        Deltamin=np.linalg.det(gammaAmin)+np.linalg.det(gammaBmin)+2*np.linalg.det(gammaCmin)
        mu1min=(.5*(Deltamin+(Deltamin**2-4*np.linalg.det(gammamin))**.5))**.5
        mu2min=(.5*(Deltamin-(Deltamin**2-4*np.linalg.det(gammamin))**.5))**.5
        kmin=np.linalg.eig(gammaAB_Bmin)
        mu3min=k[0]
        mu3min=mu3min[0]
        n1min=(mu1min-1)/2
        n2min=(mu2min-1)/2
        n3min=(mu3min-1)/2
        Smin[i]=g(n1min)+g(n2min)-g(n3min)
    
    STRING='$\epsilon$=%s' % xi
    K=.8*Iab-S
    Kmin=.8*Iabmin-Smin
    plt.semilogy(d,K,label=STRING)
    plt.semilogy(d,Kmin,color='red')


#Lower10=np.array(KS10-KS10min)
#Upper10=np.array(KS10max-KS10)
plt.semilogy(DS10,KS10,marker='.',linestyle='None',markersize=10,color='blue')
plt.semilogy(DS11,KS11,marker='.',linestyle='None',markersize=10,color='green')
plt.axis([0,30,1e-6,1e0])
plt.xlabel('Transmission distance [km]')
plt.ylabel('Secret key rate')
plt.legend(loc=1)
#plt.show()
plt.savefig('keyrate.pdf')