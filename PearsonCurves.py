import numpy as np
import scipy.stats
import scipy.special
import scipy.integrate
from scipy.integrate import quad
import pandas as pd
import math
import matplotlib.pyplot as plt
mycsv=pd.read_csv('data1.csv')
c=mycsv.X
x=c
print "Do you want to divide by class height? If yes, Enter class height, else, enter 0"
ch=input()
if ch:
    #We divide by the class height(basically) so that Xi+1-Xi=1
    x = [d / ch for d in c]
n=len(x)
print "Number of data points = ",n
#-----------------------------------------------------------------------
b=mycsv.Y
N=np.sum(b)
f=b
#sum(f)=1, makes calculation of mean simpler
#f = [e / N for e in b] 

#CALCULATING MEAN------------------------------------------------------
mean=0.0 
i=0 
for i in range(n):
    mean=mean+(f[i]*x[i])
mean=mean/N
#CALCULATING MODE------------------------------------------------------
i=0
for i in range(n):
    if(np.max(f)==f[i]):
        break;
mode=x[i]
print "Mean = ",mean
print "Mode = ",mode

#Calculating Moments about mean-----------------------------------------
mom1 = [(z-mean)**(1) for z in x]
mom2 = [(z-mean)**(2) for z in x]
mom3 = [(z-mean)**(3) for z in x]
mom4 = [(z-mean)**(4) for z in x]

mom1= np.sum(mom1*np.transpose(f))/N
mom2= np.sum(mom2*np.transpose(f))/N
mom3= np.sum(mom3*np.transpose(f))/N
mom4= np.sum(mom4*np.transpose(f))/N
print "1st moment ",mom1,"\n2nd moment ",mom2,"\n3rd moment ",mom3,"\n4th moment ",mom4
beta1=(mom3**(2))/(mom2**(3))
beta2=mom4/(mom2**2)
print "Beta1 = ",beta1,"\nBeta2 = ",beta2
kappa=(beta1*((beta2+3)**2))/((4*beta2-3*beta1)*(2*beta2-3*beta1-6)*4)
print kappa
#-----------------------------------------------------------------CHECK
arb1=5*beta2-6*beta1-9
arb2=2*beta2-3*beta1-6
Ptype=[]
for i in range(12):
    Ptype.append(0)
if kappa<0:
    print "Type 1"
    Ptype[0]=1
    if kappa>-0.15:#kappa~0
        if (beta1>-0.15)&(beta1<0.15)&(beta2>2.85)&(beta2<3.15):#beta1~0 , beta2~3
            print "Normal"
            Ptype[1]=1
        elif (beta1>-0.15)&(beta1<0.15)&(beta2<=2.85):#beta1~0 , beta2<3
            print "Type 2"
            Ptype[1]=1
        elif (beta1>-0.15)&(beta1<0.15)&(beta2>=3.15):#beta1~0 , beta2>3
            print "Type 7"
            Ptype[6]=1
        else:
            if arb1<0:
                print "type 8"
                Ptype[7]=1
            elif (arb1>0)&(arb2<0):
                print "type 9"
                Ptype[8]=1
            else:
                Ptype[0]=1
    else:
        if abs(arb2)<0.999:#doubt, kappa~INF
            print "type 3"
            Ptype[2]=1
        else:
            if arb1<0:
                print "type 8"
                Ptype[7]=1
            elif (arb1>0)&(arb2<0):
                print "type 9"
                Ptype[8]=1
            else:
                Ptype[0]=1
elif (kappa>0)&(kappa<1):
    print "type 4"
    Ptype[3]=1
    if kappa<0.15:
        if (beta1>-0.15)&(beta1<0.15)&(beta2>2.85)&(beta2<3.15):
            print "Normal"
            Ptype[1]=1
        elif (beta1>-0.15)&(beta1<0.15)&(beta2<=2.85):
            print "Type 2"
            Ptype[1]=1
        elif (beta1>-0.15)&(beta1<0.15)&(beta2>=3.15):
            print "Type 7"
            Ptype[6]=1
        else:
            Ptype[3]=1
    elif kappa>0.74:#kappa~1
        if arb2>0:
            print "type 11 (or 5?)"#ISSUE
            Ptype[10]=1
            Ptype[4]=1
        else:
            print "type 5"
            Ptype[4]=1
    else:
        Ptype[3]=1
else:
    Ptype[5]=1
    print "type 6"
    if kappa<1.15:
        print "type 5"
        Ptype[4]=1
    if arb2>0:
        print "type 11"            
        Ptype[10]=1
    else:
        Ptype[5]=1     
#-----------------------------------------------------------Finding Constants
if(Ptype[0]):
#TYPE 1
    r=6*(beta2-beta1-1)/(3*beta1-2*beta2+6)
    b=0.5*(mom2*(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)
    if mom3>0:
        m2=0.5*(r-2+(r*(r+2))*((beta1/(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)))
        m1=0.5*(r-2-(r*(r+2))*((beta1/(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)))
    else:
        m1=0.5*(r-2+(r*(r+2))*((beta1/(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)))
        m2=0.5*(r-2-(r*(r+2))*((beta1/(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)))
    a2=m2*b/(m1+m2)
    a1=m1*b/(m1+m2)
    #print "r=",r,",b=",b,",m1=",m1,",m2=",m2,",a1=",a1,",a2=",a2
    if (m1>0)&(m2>0):
        y0=N*((m1**m1)*(m2**m2)*scipy.special.gamma(m1+m2+2))/(b*((m1+m2)**(m1+m2))*scipy.special.gamma(m1+1)*scipy.special.gamma(m2+1))
        plx = np.linspace(-10.0, 100.0,100000)
        ply = y0*((1+((plx-mode)/a1))**(m1))*((1-((plx-mode)/a2))**(m2))
        #x2 = [z-mode for z in x]
        x2=x
        plt.plot(plx,ply,'r',label='Pearson\'s Curve')
    else:#changing Origin to start of curve otherwise there will be issues with 'powering'
        y0=N*(scipy.special.gamma(m1+m2+2))/((b**(m1+m2+1))*scipy.special.gamma(m1+1)*scipy.special.gamma(m2+1))
        plx = np.linspace(-10.0, 100.0,100000)
        ply = y0*((plx-x[0])**m1)*((b-(plx-x[0]))**(m2))
        #x2 = [z-x[0] for z in x]
        x2=x
        plt.plot(plx,ply,'r')
    #print y0
#---------------------------------------------------------------------------------------------
if Ptype[1]:
#TYPE 2
    m=(5*beta2-9)/(2*(3-beta2))
    asq=(2*mom2*beta2)/(3-beta2)
    y0=(N*scipy.special.gamma(2*m+2))/((asq**(1/2.0))*(2**(2*m+1))*((scipy.special.gamma(m+1))**2))
    plx = np.linspace(-10.0, 100.0,100000)
    ply = y0*(((1-((plx-mean)**2)/asq))**m)
    #x2 = [z-mean for z in x]
    x2=x
    plt.plot(plx,ply,'g')
#---------------------------------------------------------------------------------------------
if Ptype[6]:
#TYPE 7
    m=(5*beta2-9)/(2*(-3+beta2))
    asq=(2*mom2*beta2)/(-3+beta2)
    y0=(N*scipy.special.gamma(m))/(((asq*math.pi)**(1/2.0))*(scipy.special.gamma(m-0.5)))
    plx = np.linspace(-10.0, 100.0,100000)
    ply = y0*(((1+((plx-mean)**2)/asq))**(-m))
    #x2 = [z-mean for z in x]
    x2=x
    plt.plot(plx,ply,'g')
#---------------------------------------------------------------------------------------------
if Ptype[2]:
#TYPE 3
    gamma=2*mom2/mom3
    a=(2*(mom2**2)/mom3)-(mom3/(2*mom2))
    p=gamma*a
    if p>0:
        y0=N*(p**(p+1))/(a*math.exp(p)*scipy.special.gamma(p+1))
        plx = np.linspace(-10.0, 20.0,100000)
        ply = y0*((1+((plx-mode)/a))**(p))*math.exp(-gamma*(plx-mode))
        #x2 = [z-mode for z in x]
        x2=x
    else:
        y0=N*gamma*((p+1)**p)/(math.exp(p+1)*scipy.special.gamma(p+1))
        a=(p+1)/gamma
        plx = np.linspace(-10, 20.0,100000)
        plx1=np.exp(-gamma*(plx-mean))
        ply = y0*((1+((plx-mean)/a))**(p))*plx1
        #x2 = [z-mean for z in x]
        x2=x
    plt.plot(plx,ply,'b')
#---------------------------------------------------------------------------------------------
if (Ptype[3]):
#TYPE 4
    r=6*(beta2-beta1-1)/(-3*beta1+2*beta2-6)
    m=(r+2)/2.0
    nu=-(r*(r-2))*((beta1/(-(beta1*((r-2)**2))+16*(r-1)))**(1/2.0))
    if mom3>0:
        if nu>0:
            nu=-nu
    else:
        if nu<0:
            nu=-nu
    print nu
    a=((-(beta1*((r-2)**2))+16*(r-1))**(1/2.0))*(mom2**(1/2.0))/4.0
    fun=lambda d:((np.sin(d))**r)*(np.exp(nu*d))
    y0m1=scipy.integrate.quad(fun,0.0,math.pi)
    y0m2=np.exp(-nu*math.pi/2.0)*a
    y0=N/(y0m1[0]*y0m2)
    plx = np.linspace(-15, 30,100000)
    ply = y0*((1+((plx-(mean+(nu*a/r)))*(plx-(mean+(nu*a/r)))/(a*a)))**(-m))*np.exp(-nu*np.arctan((plx-(mean+(nu*a/r)))/a))
    #x2 = [z-(mean+(nu*a/r)) for z in x]
    x2=x
    plt.plot(plx,ply,'b')
#---------------------------------------------------------------------------------------------
if (Ptype[4]):
#TYPE 5
    p= 4+((8+4*(4+beta1)**(1/2.0))/(beta1))
    gamma=(p-2)*((mom2*(p-3))**(1/2.0))
    print mom3,gamma, p
    if mom3>0:
        if gamma<0:
            gamma=-gamma
    else:
        if gamma>0:
            gamma=-gamma
    print gamma
    if gamma>0:
        y0=N*(gamma**(p-1))/scipy.special.gamma(p-1)
        plx = np.linspace(-15, 30,100000)
        ply=y0*(plx-x[0])**(-p)*np.exp(-gamma/(plx-x[0]))
        #x2=[z-x[0] for z in x]
    else:
        y0=N*((p-2)**p)/(gamma*scipy.special.gamma(p-1)*np.exp(p-2))
        A=gamma/(p-2)
        plx = np.linspace(-30, 30,100000)
        ply=y0*((1+((plx-mean)/A))**(-p))*np.exp((p-2)/(1+((plx-mean)/A)))
        print y0,A,p
    x2=x    
    plt.plot(plx,ply,'r')
#---------------------------------------------------------------------------------------------
if (Ptype[5]):
#TYPE 6
    r=6*(beta2-beta1-1)/(3*beta1-2*beta2+6)
    a=0.5*(mom2*(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)
    if mom3<0:
        if a>0:
            a=-a
    else:
        if a<0:
            a=-a
    q1=-0.5*(r-2-(r*(r+2))*((beta1/(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)))
    q2=0.5*(r-2+(r*(r+2))*((beta1/(beta1*((r+2)**2)+16*(r+1)))**(1/2.0)))
    if q1<q2:
        temp=q1
        q1=q2
        q2=temp
    y0=N*(a**(q1-q2-1))*scipy.special.gamma(q1)/(scipy.special.gamma(q2+1)*scipy.special.gamma(q1-q2-1))
    plx = np.linspace(-15, 15,100000)
    ply=y0*(((plx-(mean-(a*(q1-1))/(q1-q2-2)))-a)**(q2))*(((plx-(mean-(a*(q1-1))/(q1-q2-2))))**(-q1))
    #x2 = [z-(mean-(a*(q1-1))/(q1-q2-2)) for z in x]
    x2=x
    plt.plot(plx,ply,'k')
#---------------------------------------------------------------------------------------------
if (Ptype[7]):
#TYPE 8
    p=[4-beta1,9*beta1-12,-24*beta1,16*beta1]
    roots=np.roots(p)
    arr=[r for r in roots if (r>=0)&(r<=1)]
    m=arr[0]
    a=(2-m)*((mom2*(3-m)/(1-m))**(1/2.0))
    if mom3>0:
        if a>0:
            a=-a
    else:
        if a<0:
            a=-a
    y0=N*(1-m)/abs(a)
    print m,a,y0
    plx = np.linspace(0, 5.7,100000)
    ply=y0*((1+((plx)/a))**(-m))
    #x2 = [z-x[n-1] for z in x]
    x2=x
    plt.plot(plx,ply,'k')
#---------------------------------------------------------------------------------------------
if (Ptype[8]):
#TYPE 9
    m_=2*arb1/arb2
    m_=m_*(4-beta1)+9*beta1-12
    m=-(24*beta1-((24*beta1)**2-64*beta1*m_)**(1/2.0))/(2*m_)
    print m
    a=(mom2**(1/2.0))*(2+m)*(((3+m)/(1+m))**(1/2.0))
    y0=N*(1+m)/a
    print m,a,y0
    plx = np.linspace(-50, 0,100000)
    ply=y0*((1+((plx-x[n-1])/a))**(m))
    #x2 = [z-x[n-1] for z in x]
    x2=x
    plt.plot(plx,ply,'b')
#---------------------------------------------------------------------------------------------
plt.plot(x2,f,'o', label='sample observations')
plt.title("Pearson's Frequency Curves")
plt.show()