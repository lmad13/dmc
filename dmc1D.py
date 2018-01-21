import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
class wavefunction:
    
    def __init__(self,nWalkers,potential,plotting=False):
        print "initialized",nWalkers," coordinates for",potential
        self.potential=potential
        self.xcoords=np.zeros((nWalkers))

        self.dtau=self.set_dtau()
        self.D=.5
        self.mass=self.set_mass('H2')
        print 'reduced mass is', self.mass
        self.sigma_dx=(2*self.D*self.dtau/self.mass)**0.5
        self.mu_dx=0.0
        self.alpha=.25/self.dtau
        self.plotting=plotting
        self.recrossing=True if potential=='half harmonic' else False

    def set_dtau(self):
        print 'set dtau to be ',5.0
        return 5.0

    def set_mass(self,molecule):
        conversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839    in Atomic Units!!
        massH=1.00782503223
        if molecule=='H2':

            return (massH*massH)/(massH+massH)*conversionFactor

    def V(self,x):
        #input in bohr, output in a.u.
        self.omega=2000.0
        k=(self.omega*2.0*np.pi*3.0*10**(10))**2*self.mass #cm-1
        convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
        k=k*convfactor
        if self.potential=='harmonic':
            v=0.5*k*(x*x)
        if self.potential=='half harmonic':
            #v=0.5*k*(x*x) if x>0 else 100000
            inf=10000.0
            v=0.5*k*(x*x)
            mask=(x<0.0)
            v[mask]=inf
            #v=[vb if c else vinf for (c,vb,vinf) in zip(x>0.0,0.5*k*(x*x),inf*np.ones(x.size))]
        return v


#function propagate
    def propagate(self,x,nSteps):

        N_size_step=x.size
        nSize=x.size
        v_ref=np.average(self.V(x))
        vRefList=[]
        
        #print 'propagating!'
        for step in range(nSteps):
            dx=self.diffuse(x)
            x=x+dx
            v=self.V(x)
            if self.plotting and step%(nSteps/10)==0:
                plt.scatter(x,v)
                xrange= np.linspace(np.min(x), np.max(x), num=5)
                plt.plot(xrange,np.zeros((5))+v_ref)
                plt.show()
            #Elimintation of a random selection of walkers in
            #the classically forbidden region
            N_r=np.random.random(N_size_step)
            P_d=1-np.exp(-(v-v_ref)*self.dtau)
            Diff=N_r-P_d
            mask_survive = (Diff>0)
            nDeaths=np.sum(np.array(Diff<0).astype(int))
            survivors=x[mask_survive]
            #print survivors.shape, 'before and after',
            #Recrossing correction will go here
            #print'Census: Deaths:',nDeaths,
            if self.recrossing:
                #print 'recrossing correction activated!!'
                P_recrossDeath=np.exp(-2.0*(x-dx)*x*np.sqrt(self.mass*self.mass)/self.dtau) ##mass is reduced mass!
                #print P_recrossDeath[0:10], x[0:10],(x-dx)[0:10], np.sqrt(self.mass*self.mass)/self.dtau
            
                if self.plotting:
                    plt.scatter(x,P_recrossDeath)
                    plt.show()
                #print 'P_recrossDeath 1', P_recrossDeath
                #P_recrossDeath=P_recrossDeath[mask_survive]
                Diff=N_r-P_recrossDeath
                #Diff=N_r[mask_survive]-P_recrossDeath  ##Should this be a different random number?
                mask_survive_recross=(Diff>0)
                tempRecrossCensus=np.sum(np.array(Diff<0).astype(int))
                #print mask_survive[0:10], mask_survive_recross[0:10]
                mask_survive=np.logical_and(mask_survive, mask_survive_recross)
                survivors=x[mask_survive]            
                #print mask_survive_recross[0:10]
                #print survivors.shape
                #print 'Recross Deaths', tempRecrossCensus,
            #Creation of a random selection of walkers in the
            P_exp_b=np.exp(-(v-v_ref)*self.dtau)-1.0
            
            weight_P_b=P_exp_b.astype(int)
            P_b=P_exp_b-weight_P_b            #classically allowed region
            #P_b=-(v-v_ref)*self.dtau

            Diff=N_r-P_b
            mask_b = (Diff<0)
            next_gen=x[mask_b] 

            nBirths=np.sum(np.array(Diff<0).astype(int))
            addBirthtot=0
            new_pop=next_gen
            for n,(particle,weight) in enumerate(zip(x,weight_P_b)):
                if weight>0: #i.e. the dead can't reproduce
                    #print 'weight:',weight
                    if weight>10:
                        #print 'weight is too big, resetting to 10'
                        #print x[n],v(n),'<',v_ref, -(v(n)-v_ref)
                        weight=10
                    addBirthtot=addBirthtot+weight

                    temp=np.tile(particle,weight)
                    new_pop=np.concatenate((new_pop,temp))
            next_gen=new_pop


            #collect survivors and next generation
            new_population=np.concatenate((survivors,next_gen))
            N_size_step=new_population.size
            #print '. Births:',nBirths, ". Add' births: ", addBirthtot
            #readjust V_ref
            v_average=np.average(self.V(new_population))
            v_ref=v_average+(self.alpha*(1-float(N_size_step)/float(nSize)))
            #print '(',N_size_step,'/',nSize,') v_ref',v_ref, '=', v_average,'+', (self.alpha*(1-float(N_size_step)/float(nSize)))
            #print 'size',N_size_step, 'compared to ', nSize
            if v_ref<0:
                print' step:',step, v_ref, ':',float(N_size_step)/float(nSize), '=', float(N_size_step),'/',float(nSize)
            vRefList.append(v_ref)
            

            x=new_population
#propagate wavefunction for nsteps WITHOUT kmnowing the final v_ref
#take in walker positions, V, nsteps, delta t
#return vref measurment (avearged over 1000 steps, only at equilibrium), snapshots of wavefunction every 500 t steps after equilibration
        print 'Average from ',nSteps/2,' steps until the end is ',np.average(vRefList[nSteps/2:])*au2wn,
        print 'final number of ancestors', N_size_step
        return np.average(vRefList[nSteps/2:])*au2wn

#function diffuse
    def diffuse(self,x):
        N_size_step=x.shape[0]
        dx=np.random.normal(self.mu_dx, self.sigma_dx, N_size_step)
        return dx

#function CalcDesc
#calculate descendant weights
#take in walker positions, V, nsteps, delta t         
#return the number of descendants



