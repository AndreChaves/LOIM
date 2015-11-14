import numpy as np
import scipy as sp
import susc_gg as sc
import scipy.optimize as opt

def test(q,f,system_parameters=[]):
    return f+q    

class parameters():
    def __init__(self,t1,t2=[]):
          self.tau=t1	
	

class Susceptibility:
    def __init__(self,formula,freq_vec,q_vec,system_parameters=[]):
            n_freq=len(freq_vec)
            n_q=len(q_vec)
            self.system_parameters=system_parameters
            self.n_freq=n_freq
            self.n_q=n_q
            self.freq=freq_vec
            self.q=q_vec
            self.formula=formula
            self.susc=np.ndarray(shape=(n_freq,n_q),dtype=complex)
            for j in range(n_freq):
                for j2 in range(n_q):
                    freq=freq_vec[j]
                    q=q_vec[j2]
                    self.susc[j,j2]=formula(q,freq,system_parameters)
            
    def Mermim_Formula(self,qin,fin,jin=[]):
        if self.system_parameters:
            g1=self.system_parameters.gamma
            if(g1!=0):
                t1=1./g1
                if(jin): #vector input
                    chi=self.susc[:,jin]
                    q=self.q[jin]
                    f=self.freq
                else: #function
                    chi=self.formula(qin,fin,self.system_parameters)
                    f=fin
                    q=qin    
                chi0=Susceptibility(self.formula,freq_vec=[0.],q_vec=[q],system_parameters=self.system_parameters)
                num=1.-1.j*t1*chi0.susc[0,0]*np.multiply(f,chi)
                den=-1.j*t1*f*chi0.susc[0,0]+chi
                return np.divide(num,den)
            else:
                return self.formula(qin,fin)
        else:
            return self.formula(qin,fin)
        
                   
    def New_Mermim_Formula(self,qin,fin,jin=[]):
        if self.system_parameters:
            g1=self.system_parameters.gamma
            if(g1!=0):
                t1=1./g1
                if(jin): #vector input
                    chi=self.susc[:,jin]
                    q=self.q[jin]
                    f=self.freq
                else: #function
                    chi=self.formula(qin,fin,self.system_parameters)
                    f=fin
                    q=qin    
                #chi0=Susceptibility(self.formula,freq_vec=[0.],q_vec=[q],system_parameters=self.system_parameters)
                #num=1.-1.j*t1*chi0.susc[0,0]*np.multiply(f,chi)
                #den=-1.j*t1*f*chi0.susc[0,0]+chi
                return np.divide(num,den)
            else:
                return self.formula(qin,fin)
        else:
            return self.formula(qin,fin) 
               
    def mtest(self):        
        def method(qin,fin,parameters,chi=[]):
            t1=1./(parameters.gamma)
            chi0=Susceptibility(self.formula,freq_vec=[0.],q_vec=[qin],system_parameters=self.system_parameters)
            if not(chi):
                    chi=self.formula(qin,fin,self.system_parameters)            
            num=chi0.susc*np.multiply(1.-1.j*t1*fin,chi)
            
            den=-1.j*t1*fin*chi0.susc+chi
            return np.divide(num,den)[0][0]
        if self.system_parameters:
            g1=self.system_parameters.gamma
            if(g1!=0):
                result=method
            else:
                result=self.formula
        else:
            result=self.formula
        return result
    
    def Mermim(self):
        try: 
                if self.system_parameters:                    
                    g1=self.system_parameters.gamma
                    if(g1!=0.):                        
                        t1=1./g1
                        freq_vec=self.freq
                        q_vec=self.q
                        freq_null=np.zeros(self.n_freq,dtype=complex)
                        self.system_parameters.gamma=0.
                        static=Susceptibility(self.formula,freq_null,q_vec,self.system_parameters)
                        self.system_parameters.gamma=g1
                        y=np.ndarray(shape=(self.n_freq,self.n_q),dtype=complex)
                        num=np.ndarray(shape=(self.n_freq,self.n_q),dtype=complex)
                        den=np.ndarray(shape=(self.n_freq,self.n_q),dtype=complex) 
                        for j in range(self.n_q):
#                            y[:,j]=Susceptibility.Mermim_Formula(self,jin=j)
                            num[:,j]=np.multiply(1.-1.j*freq_vec*t1,np.multiply(static.susc[:,j],self.susc[:,j]))
                            den[:,j]=-1.j*np.multiply(freq_vec*t1,static.susc[:,j])+self.susc[:,j]	
                        y=np.divide(num,den)
                    else:
                        print("Warning: tau=infinity")
                        y=self.susc                        
                else:
                    print("Warning: tau=infinity")
                    y=self.susc   
        except ValueError:
                print ("Warning: some problem ocurred")  
                y=self.susc
        self.merm=y
#        return y
    def Conductivity(self,Mermim_option=False):         
        y=np.ndarray(shape=(self.n_freq,self.n_q),dtype=complex)
        if(Mermim_option):
            print(self.system_parameters)
            self.Mermim()
            z=self.merm
        else:
            z=self.susc
        for j in range(self.n_q):
            if(self.q[j]!=0.):
                y[:,j]=1.j*np.divide(np.multiply(self.freq,z[:,j]),np.power(self.q[j],2.))
        return y
    def Plasmon(self,q_vec,Potential=[],step=1.E-3,xend=2.,Mermim_option=False):
        system_parameters=self.system_parameters
        if not(Potential):
            Potential= lambda x: 1./(x*x)
        if not(Mermim_option):                
            susc_formula=self.formula             
        else:
            susc_formula=self.mtest()
#        q=self.q[0]
        n_q=len(q_vec)        
        xplot=np.ndarray(shape=(n_q),dtype=complex)
        xtry=np.ndarray(shape=(n_q),dtype=float)
#        f_guess=q_vec[0]*1.1
####First Zero        
        zeros=Susceptibility.first_zero( Susceptibility.P_zero_real,x_start=q_vec[0]+1.E-9,step_size=step,x_end=xend,args=(q_vec[0],Potential,susc_formula,system_parameters))
######
        xvec=[]
        xvec2=[]
        for n_zero in range(len(zeros)):
            f_guess=zeros[n_zero]
            f_m=0.5*(f_guess[0]+f_guess[1])
            for j in range(n_q-1):
                q=q_vec[j+1]                
                xreal=opt.newton(Susceptibility.P_zero_real,f_m,args=(q,Potential,susc_formula,system_parameters))
                f_m=xreal
                xtry[j]=xreal
#            print(xreal,Susceptibility.P_zero_real(xreal,q,Potential,susc_formula,system_parameters))
                qi=[xreal,0.]#np.array([1.,0.])
                xout=opt.root(Susceptibility.P_zero,qi,args=(q,Potential,susc_formula,system_parameters),method='lm')
#            print(xout.x)
                xplot[j]=xout.x[0]+1.j*xout.x[1]
#            xplot[j]=xreal
            xvec.append(xtry)
            xvec2.append(xplot)
        return xvec,xvec2#xtry#xplot,xtry
#        print(xout,Susceptibility.P_zero_real([xreal,q,Potential,susc_formula,system_parameters)))
#        vetor_provisorio=np.ndarray(shape=(self.n_freq),dtype=float)
#        for j in range(self.n_freq):
#            vetor_provisorio[j]=Susceptibility.P_zero_real(self.freq[j],q,Potential,susc_formula,system_parameters)
#        return vetor_provisorio
    def P_zero(ftotal,q,Potential,susc_formula,system_parameters): 
            f=ftotal[0]+1.j*ftotal[1]
            valor=1.-Potential(q)*susc_formula(q,f,system_parameters)      
            return valor.real,valor.imag   
        
    def P_zero_real(f,q,Potential,susc_formula,system_parameters):    
            valor=1.-Potential(q)*susc_formula(q,f,system_parameters)      
            return valor.real              

    def first_zero( func,x_start,x_end,step_size=1.E-3,args=()):        
            n_interact=int((x_end-x_start)/step_size)
            nzeros=0
            x0=x_start
            x1=x0+step_size
            f0=func(*((x0,) + args))
            zerop=[]
            for j in range(n_interact):
                f1=func(*((x1,) + args))
                if(f0*f1<0):
                        nzeros=nzeros+1
                        zerop.append([x0,x1]) 
                x0=x1
                f0=f1
                x1=x1+step_size
            if(nzeros==0):
                print("warning: there's no zeros")
            zeros=np.asarray(zerop)
            return zeros 
#class unit_convert():
#    def ev2hz():
#        
#    def ev


