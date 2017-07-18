from django.shortcuts import render
from django.http import HttpResponse
import numpy as np
import pandas as pd
import math
import scipy.stats
import scipy.special
import scipy.integrate
import matplotlib.pyplot as plt
import io,StringIO
import base64
def index(request):
    context = {}
    return render(request, 'polls/index.html', context)

def result(request):
    if request.POST and request.FILES:
        mycsv=pd.read_csv(request.FILES['csv_file'])
        c=mycsv.X
        #We divide by the class height(basically) so that Xi+1-Xi=1
        x = [d / 5 for d in c]
        n=len(x)
        b=mycsv.Y
        N=np.sum(b)
        f = b 
        #Thus sum(f)=1, makes calculation of mean simpler

        #CALCULATING MEAN
        mean=0.0 
        i=0 
        for i in range(n):
            mean=mean+f[i]*x[i] 
        mean=mean/N

        #CALCULATING MODE
        i=0
        for i in range(n):
            if(np.max(f)==f[i]):
                break
        mode=x[i]

        #Calculating Moments about mean
        mom1 = [(z-mean)**(1) for z in x]
        mom2 = [(z-mean)**(2) for z in x]
        mom3 = [(z-mean)**(3) for z in x]
        mom4 = [(z-mean)**(4) for z in x]

        mom1= np.sum(mom1*np.transpose(f))/N
        mom2= np.sum(mom2*np.transpose(f))/N
        mom3= np.sum(mom3*np.transpose(f))/N
        mom4= np.sum(mom4*np.transpose(f))/N

        #Calculating beta1 and beta2
        beta1=(mom3**(2))/(mom2**(3))
        beta2=mom4/(mom2**2)
        #Calculating kappa criterion
        kappa=(beta1*((beta2+3)**2))/((4*beta2-3*beta1)*(2*beta2-3*beta1-6)*4)

        #Performing Checks
        arb1=5*beta2-6*beta1-9
        arb2=2*beta2-3*beta1-6
        Ptype=[]
        for i in range(12):
            Ptype.append(0)
        if kappa<0:
            Ptype[0]=1
            if kappa>-0.15:#kappa~0
                if (beta1>-0.15)&(beta1<0.15)&(beta2>2.85)&(beta2<3.15):#beta1~0 , beta2~3
                    Ptype[1]=1
                elif (beta1>-0.15)&(beta1<0.15)&(beta2<=2.85):#beta1~0 , beta2<3
                    Ptype[1]=1
                elif (beta1>-0.15)&(beta1<0.15)&(beta2>=3.15):#beta1~0 , beta2>3
                    Ptype[6]=1
                else:
                    if arb1<0:
                        Ptype[7]=1
                    elif (arb1>0)&(arb2<0):
                        Ptype[8]=1
                    else:
                        Ptype[0]=1
            else:
                if abs(arb2)<0.999:#doubt, kappa~INF
                    Ptype[2]=1
                else:
                    if arb1<0:
                        Ptype[7]=1
                    elif (arb1>0)&(arb2<0):
                        Ptype[8]=1
                    else:
                        Ptype[0]=1
        elif (kappa>0)&(kappa<1):
            Ptype[3]=1
            if kappa<0.15:
                if (beta1>-0.15)&(beta1<0.15)&(beta2>2.85)&(beta2<3.15):
                    Ptype[1]=1
                elif (beta1>-0.15)&(beta1<0.15)&(beta2<=2.85):
                    Ptype[1]=1
                elif (beta1>-0.15)&(beta1<0.15)&(beta2>=3.15):
                    Ptype[6]=1
                else:
                    Ptype[3]=1
            elif kappa>0.74:#kappa~1
                if arb2>0:#doubt
                    Ptype[10]=1
                    Ptype[4]=1
                else:
                    Ptype[4]=1
            else:
                Ptype[3]=1
        else:
            Ptype[5]=1
            if kappa<1.15:
                Ptype[4]=1
            if arb2>0:           
                Ptype[10]=1
            else:
                Ptype[5]=1     
        
        if(Ptype[0]):
        #TYPE 1
        #Huge Problem of values going to negetive
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
                plt.plot(plx,ply,'r')
            else:#changing Origin to start of curve otherwise there will be issues with 'power'
                y0=N*(scipy.special.gamma(m1+m2+2))/((b**(m1+m2+1))*scipy.special.gamma(m1+1)*scipy.special.gamma(m2+1))
                plx = np.linspace(-10.0, 100.0,100000)
                ply = y0*((plx-x[0])**m1)*((b-(plx-x[0]))**(m2))
                #x2 = [z-x[0] for z in x]
                x2=x
                plt.plot(plx,ply,'r')
                #print y0

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
        

        #PLOTTING GRAPH
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(plx,ply)
        ax.plot(x2, f,'o')
        plt.title('Pearson\'s Curves')
        img_in_memory = StringIO.StringIO()
        fig.savefig(img_in_memory, format="png")
        image=base64.b64encode(img_in_memory.getvalue()) 

        context={'image':image,'mean':mean,'mode':mode,'beta1':beta1,'beta2':beta2,'kappa':kappa,}
        img_in_memory.close()
        return render(request, 'polls/result.html', context)
    else:
        return HttpResponse("Form Not Submitted")