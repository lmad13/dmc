import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
class wavefunction:
    
    def __init__(self,nWalkers,potential,plotting=False,omegaInput=2000.00,molecule='H2'):
        print "initialized",nWalkers," coordinates for",potential
        self.potential=potential
        self.xcoords=np.zeros((nWalkers))
        self.dtau=self.set_dtau()
        self.D=.5
        self.mass=self.set_mass(molecule)
        self.omega=omegaInput
        print 'reduced mass is', self.mass
        self.sigma_dx=(2.0*self.D*self.dtau/self.mass)**0.5
        print 'sigma_dx',self.sigma_dx
        self.mu_dx=0.0
        self.alpha=0.500/self.dtau
        self.plotting=plotting
        self.recrossing=True if 'half harmonic' in potential else False
        print 'recrossing is ', self.recrossing
    def setX(self,x):
        self.xcoords=x

    def set_dtau(self):
        print 'set dtau to be ',10.0
        return 10.0
    def get_dtau(self):
        
        return self.dtau

    def set_mass(self,molecule):
        #why? for reduced mass in amu
        conversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839    in Atomic Units!!
        massH=1.00782503223
        massLi=7.0
        if molecule=='H2':
            massAtom=massH
        elif molecule=='Li2':
            massAtom=massLi
        return (massAtom*massAtom)/(massAtom+massAtom)*conversionFactor

    def getIdealVref(self):
        if self.potential=='harmonic':
            return self.omega/(2.0*au2wn)
        
        elif 'half harmonic' in self.potential:
            return 3.0*self.omega/(2.0*au2wn)
        
    def getAlpha(self):
        k=(self.omega*2.0*np.pi*3.0*10**(10))**2*self.mass #cm-1                                            
        convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2   
        k=k*convfactor
        #print 'my k is ', k
        alpha=np.sqrt(k*self.mass)
        return alpha

    def V(self,x):
        #input in bohr, output in a.u.
        k=(self.omega*2.0*np.pi*3.0*10**(10))**2*self.mass #cm-1
        convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
        k=k*convfactor
        if self.potential=='harmonic':
            
            v=0.5*k*(x*x)
        if self.potential=='half harmonic right':
            #v=0.5*k*(x*x) if x>0 else 100000
            inf=10000.0
            v=0.5*k*(x*x)
            mask=(x<0.0)
            v[mask]=inf
            #v=[vb if c else vinf for (c,vb,vinf) in zip(x>0.0,0.5*k*(x*x),inf*np.ones(x.size))]
        if self.potential=='half harmonic left':
            #v=0.5*k*(x*x) if x>0 else 100000
            inf=10000.0
            v=0.5*k*(x*x)
            mask=(x>0.0)
            v[mask]=inf
        elif self.potential=='half harmonic':
            v=0.5*k*(x*x)

        return v


#function propagate
    def propagate(self,x,nSteps,setV_ref=False,ConstantV_ref=0,printCensus=False,nSize=0,plotWalkers=False):

        descendants=np.zeros((x.size))  #not yet implemented
        whoYaFrom=np.arange(x.size)
        N_size_step=x.size
        if nSize==0:#if it is the default value of zero...it needs to be set.
            nSize=x.size
        if setV_ref:
            v_ref=ConstantV_ref
        else:
            v_ref=np.average(self.V(x))#+(self.alpha*(1-float(N_size_step)/float(nSize)))
        vRefList=[]
        population=[]
        if plotWalkers:
            plt.figure(2)

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
            
            if printCensus: print'Census: Deaths:',nDeaths, 

            #recrossing correction for pop near node, only true for half harmonic PES
            if self.recrossing:
                crossed=((x-dx)*x<0)
                P_recrossDeath=np.exp(-2.0*(x-dx)*x*np.sqrt(self.mass*self.mass)/self.dtau) ##mass is reduced mass!
                if self.plotting:
                    plt.scatter(x,P_recrossDeath)
                    plt.show()
                N_r_recross=np.random.random(N_size_step)
                Diff=N_r_recross-P_recrossDeath
                mask_survive_recross=(Diff>0.0000000000)
                tempRecrossCensus=np.sum(np.array(Diff<0).astype(int))
                mask_survive=np.logical_and(mask_survive, mask_survive_recross)
                survivors=x[mask_survive]            

            #Creation of a random selection of walkers in the classically allowed region
            P_exp_b=np.exp(-(v-v_ref)*self.dtau)-1.0
            P_exp_b[np.logical_not(mask_survive)]=0.000000 #BECAUSE THE DEAD CANNOT REPRODUCE!! NO ZOMBIE MOTHERS!! 
            if printCensus: print '\n P_exp_b',np.average(P_exp_b[(P_exp_b>0)]),np.std(P_exp_b[(P_exp_b>0)])
            weight_P_b=P_exp_b.astype(int)
            P_b=P_exp_b-weight_P_b            #classically allowed region
            #P_b[np.logical_not(mask_survive)]=0.0
            Diff=N_r-P_b
            mask_b = (Diff<0)
            next_gen=x[mask_b] 
            new_pop_whoYaFrom=whoYaFrom[mask_b]
            nBirths=np.sum(np.array(Diff<0).astype(int))
            addBirthtot=0
            new_pop=next_gen
            #for the additional births
            for n,(particle,weight) in enumerate(zip(x,weight_P_b)):
                if weight>0: #i.e. the dead can't reproduce
                    if weight>10:
                        #this really shouldn't happen
                        print 'weight is too big, resetting to 10'
                        print x[n],v(n),'<',v_ref, -(v(n)-v_ref)
                        weight=10
                    addBirthtot=addBirthtot+weight

                    temp=np.tile(particle,weight)
                    temp_whoYaFrom=np.tile(whoYaFrom[p],weight)
                    new_pop=np.concatenate((new_pop,temp))
                    new_pop_whoYaFrom=np.concatenate((new_pop_whoYaFrom,temp_whoYaFrom))

            if printCensus: print '. Births:',nBirths, ". Add' births: ", addBirthtot
            if plotWalkers:
                plt.scatter(x[mask_survive],v[mask_survive],c='black')
                plt.scatter(x[mask_b],v[mask_b],c='blue',s=(P_exp_b[mask_b]*100.0)**2)
                plt.scatter(x[(np.logical_not(mask_survive))],v[(np.logical_not(mask_survive))],c='red',s=(P_d[(np.logical_not(mask_survive))]*100.0)**2)
                plt.plot([-1.0,1.0],[v_ref,v_ref],c='magenta')
                #plt.quiver(x,v,dx,np.ones(N_size_step))
                plt.show()

            #readjust V_ref
            next_gen=new_pop
            next_gen_whoYaFrom=new_pop_whoYaFrom
            #collect survivors and next generation
            new_population=np.concatenate((survivors,next_gen))
            N_size_step=new_population.size
            v_average=np.average(self.V(new_population))
            if not setV_ref:
                v_ref=v_average+(self.alpha*(1.00-float(N_size_step)/float(nSize)))
            if printCensus: print '(',N_size_step,'/',nSize,') v_ref',v_ref, '=', v_average,'+', (self.alpha*(1-float(N_size_step)/float(nSize)))

            if v_ref<0 and step>5:
                print 'this is problematic.  NSize is probably too small'
                print' step:',step, v_ref, ':',float(N_size_step)/float(nSize), '=', float(N_size_step),'/',float(nSize)

            vRefList.append(v_ref)
            population.append(N_size_step)
            whoYaFrom=np.concatenate((whoYaFrom[mask_survive],next_gen_whoYaFrom))
            x=new_population
            #print x.shape, whoYaFrom.shape
        #print 'Average from ',nSteps/2,' steps until the end is ',np.average(vRefList[nSteps/2:])*au2wn,
        #print 'final number of ancestors', N_size_step
        for anc in whoYaFrom:
            descendants[anc]=descendants[anc]+1

        return vRefList, population, x, descendants


    def diffuse(self,x):
        N_size_step=x.shape[0]
        dx=np.random.normal(self.mu_dx, self.sigma_dx, N_size_step)
        return dx

#function CalcDesc
#calculate descendant weights
#take in walker positions, V, nsteps, delta t         
#return the number of descendants



