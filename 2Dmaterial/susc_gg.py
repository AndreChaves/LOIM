#from numpy import *
import numpy as np
I=1.j
pi=np.pi


def part1(q,f):
    if(f.real>q):
        resp=-1./(4.*pi)*q*q/np.sqrt(f**2-q**2)
    else: 
        resp=-1./(4.*pi)*q*q/np.sqrt(q**2-f**2)    
    return resp

def Gless(x,x0):
    resp=x*np.sqrt(x0*x0-x*x)-(2.-x0*x0)*np.arccos(x/x0)
    return resp

def x_0(q,f,m):
    resp=np.sqrt(1+4.*m**2/(q**2-f**2))
    return resp

def x_02(q,f,m):
    resp=1+4.*m*m/(q*q-f*f)
    return resp

def Ggreater(x,x0):
    resp=x*np.sqrt(x*x-x0*x0)-(2.-x0*x0)*np.arccosh(x/x0)
    return resp

def G0(x,x0):
    resp=x*np.sqrt(x*x-x0*x0)-(2.-x0*x0)*np.arcsinh(x/np.sqrt(-x0*x0))
    return resp


def imag_pol0(q,f,m):
    resp=0.
    x=f*f-q*q-4.*m*m
    if(x.real>0):
        x02=x_02(q,f,m)
        resp=-q*q/2./np.sqrt(f*f-q*q)*(1.-0.5*x02)
    return resp

def real_pol0(q,f,m):
    resp=0.
    func=2./pi*q*q
    par=m/(q*q-f*f)
    if(q>f.real):
        z=(q*q-f*f-4.*m*m)/(f*f-q*q-4.*m*m)
        par2=np.arccos(z)
    else: 
        x=np.sqrt(f*f-q*q)
        z=(2.*m+x)/(2.*m-x)
        if(2.*m>x.real):
            par2=-np.log(z)
        else:
            par2=-np.log(-z)
    if(f.real>q):
        func2=(q*q-f*f-4.*m*m)/(4.*np.sqrt(f*f-q*q)**3)
    else:
        func2=(q*q-f*f-4.*m*m)/(4.*np.sqrt(q*q-f*f)**3)
    resp=func*(par+func2*par2)
    return resp

def imag_pol(q,f,m,mu):
    resp=0.
    x0=x_0(q,f,m)
    func=part1(q,f)
    reg=test_(q,f,m,mu)
    xm=(2.*mu+f)/q
    xl=(2.*mu-f)/q
    if(reg=='1A'):
        resp=Ggreater(xm,x0)-Ggreater(xl,x0)
    if(reg=='2A'):
        resp=Ggreater(xm,x0)
    if(reg=='2B'):
        resp=-Gless(xl,x0)
    if(reg=='3B'):
        resp=pi*(2.-x0*x0)
    if(reg=='4B'):
        resp=pi*(2.-x0*x0)
    resp=resp*func 
    return resp


def real_pol(q,f,m,mu):
    resp=0.
    x0=x_0(q,f,m)
    func=part1(q,f)
    reg=test_(q,f,m,mu)
    xm=(2.*mu+f)/q
    xl=(2.*mu-f)/q
    if(reg=='1A'):
        resp=0.
    if(reg=='2A'):
        resp=Gless(xl,x0)
    if(reg=='3A'):
        resp=Gless(xm,x0)+Gless(xl,x0)
    if(reg=='4A'):
        resp=Gless(xl,x0)-Gless(xm,x0)
    if(reg=='1B'):
        resp=Ggreater(xm,x0)-Ggreater(xl,x0)
    if(reg=='2B'):
        resp=Ggreater(xm,x0)
    if(reg=='3B'):
        resp=Ggreater(xm,x0)-Ggreater(-xl,x0)
    if(reg=='4B'):
        resp=Ggreater(xm,x0)+Ggreater(-xl,x0)
    if(reg=='5B'):
        resp=G0(xm,x0)-G0(xl,x0)
    resp=2.*mu/pi+resp*func 
    return resp

def test_(q,f_complex,m,mu):
    Ef=mu
    if(Ef>m):
        qf=np.sqrt(Ef**2-m*m)
    else:
        qf=0.
    f=f_complex.real
    reg='X'
    z1=np.sqrt((q-qf)**2+m**2)
    z2=np.sqrt((q+qf)**2+m**2)
    z3=np.sqrt(q*q+4.*m*m)
    if(q>f):
        if(f<Ef-z1):
            reg='1A'
        if(f<-Ef+z2):
            if(f>np.fabs(z1-Ef)):
                reg='2A'
        if(f<-Ef+z1):
            reg='3A'
        if(f>-Ef+z2):
            reg='4A'
    else:
        if(q<=2.*qf and f>=z3 and f<=mu+z1):
            reg='1B'
        if(f>=mu+z1 and f<=z2+mu):
            reg='2B'
        if(f>=mu+z2):
            reg='3B'
        if(q>=2.*qf and f>=z3 and f<=mu+z1):
            reg='4B'
        if(f<=z3):
            reg='5B'
    if(reg=='X'):
        print ('X',q,f)
    return reg



def test_parameters(parameters=[]):
    if parameters:
        if(hasattr(parameters,'Ef')): 
            Ef=parameters.Ef
        else:
            Ef=0.
        if(hasattr(parameters,'gamma')):   
            gamma=parameters.gamma
        else:
            gamma=0.
        if(hasattr(parameters,'lam')):   
            lam=parameters.lam
        else:
            lam=0.
        if(hasattr(parameters,'mass')):
            mass=parameters.mass
        else:
            mass=0.   
        m=mass
        delta=m
        mu=Ef
        if(hasattr(parameters,'doped')):
            doped=parameters.doped
        else:
            doped=False
    else:
        Ef=gamma=tau=lam=m=delta=mu=0.
        doped=False
    return Ef,gamma,lam,m,delta,mu,doped

def graphene_susc(q,ome,parameters=[]):
    Ef,gamma,lam,m,delta,mu,doped=test_parameters(parameters)
    y=0.   
    s=-1.
    y2=0.
    y0=0.
    y20=0.
    for j1 in range(2):
        tau=-1.
        for j2 in range(2):
            m_ef=m-lam*s*tau/2.
            mu_ef=mu+lam*s*tau/2.
            y=y+imag_pol(q,ome+1.j*gamma,m_ef,mu_ef)/4.
            y2=y2-real_pol(q,ome+1.j*gamma,m_ef,mu_ef)/4.
            y0=y0+imag_pol0(q,ome+1.j*gamma,m_ef)/4.
            y20=y20-real_pol0(q,ome+1.j*gamma,m_ef)/4.
            if(mu_ef<delta): 
                y=y0  
                y2=y20
            if(doped):
                y=y-y0
                y2=y2-y20
            tau=-tau
        s=-s
    return y2+1.j*y
 
class grap_parameters():
    def __init__(self,Ef=0.,gamma=0.,lam=0.,m=0.,doped=False):     
            self.gamma=gamma
            self.Ef=Ef
            self.lam=lam
            self.mass=m	
            self.doped=doped

def graphene_susc0(q,ome,parameters=[]):
    import graphene0 
    from numpy import sqrt
    Ef,gamma,lam,m,delta,mu,doped=test_parameters(parameters)
    if(Ef!=0.):
        x=q/Ef
        y=(ome+1.j*gamma)/Ef
        chi=Ef*graphene0.graphene_susc.polarizibility(y,x)
    else:
        q2=q*q
        f2=ome*ome
        if(q>ome):
             chi=-q2/(4.*np.sqrt(q2-f2))
        else:
             chi=-1.j*q2/(4.*np.sqrt(f2-q2))
    return chi


def graf_cond(ome,parameters=[]):
    Ef,gamma,lam,m,delta,mu,doped=test_parameters(parameters)
    drude=-np.divide(4.*Ef,1.j*ome-gamma)/pi
    inter=1.+1./pi*(np.arctan( (ome-2.*Ef)/gamma )-np.arctan((ome+2.*Ef)/gamma))
    num=np.power(2.*Ef-ome,2)+ gamma*gamma
    den=np.power(2.*Ef+ome,2)+ gamma*gamma
    inter=inter+1.j/2./pi*np.log(np.divide(num,den))
    cond=inter+drude
    return cond