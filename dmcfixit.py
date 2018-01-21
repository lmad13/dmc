print 'start dmc'
print 'YOU GOT HERE!'
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import time
#from drawnow import drawnow
#objective : Implement JBA 1975 Paper

def V_inv(fun):
    print 'V_inv omega', omega
    return fun


def initialize(potentialInfo, mu, sigma, N,omega_vibration):
    (potential,V,invPotential)=potentialInfo
    print 'initialize omega', omega_vibration
    if potential=='harmonic' or potential=='well' or potential=='box':
        a=-2.0/np.log(omega_vibration/10.0)
        #reasonableRange=(a,-a)
        reasonableRange=(-1.5,1.5)
        print 'REASONABLE RANGE:', reasonableRange
        x1=np.zeros(N)+invPotential(1.1*omega_vibration/2.0,omega_vibration)
        print 'x data!!(max, min, ave)', np.max(x1), np.min(x1), np.average(x1)
        

    elif potential=='halfHO':
#        x1=np.random.normal(mu, sigma, N)
        a=-2.0/np.log(omega_vibration/10.0)
        reasonableRange=(a,-a)
        
        print 'REASONABLE RANGE:', reasonableRange
        x1=np.zeros(N)+invPotential(1.1*omega_vibration*3.0/2.0,omega_vibration)
        print 'x data!!(max, min, ave)', np.max(x1), np.min(x1), np.average(x1)

    elif potential=='double-well' or 'halfDoubleWell':
        k=omega_vibration**2
        b=omega_vibration*2.0
        a=b/2.0
        g=b**2/(4.0*a)
        minLocations=np.sqrt(b/(2.0*a))
        xleft=np.random.normal(mu-minLocations,sigma, N/2)
        xright=np.random.normal(mu+minLocations,sigma, N/2)
        #x1=np.concatenate((xleft,xright))
        x1=np.concatenate((xleft,xleft))
        a=-1*invPotential(omega_vibration*100.0,omega_vibration) #omega_vibration*100.0,potential,omega=omega_vibration)
        reasonableRange=(a,-a)
        print 'REASONABLE RANGE:', reasonableRange
    else:
        x1=np.zeros(N)
    return x1,reasonableRange
    

def set_dtau(potential,omega_vibration,recrossing):
    print 'set_dtau based on omega', omega_vibration
    if potential=='harmonic' or potential=='halfHO':
        dtau=np.square(1/omega_vibration)  ## why did I do this?
#        if recrossing:
#            dtau=10.0*dtau
    elif potential=='double-well' or potential=='halfDoubleWell':
        dtau=np.square(1/omega_vibration)  ## why did I do this?

    elif potential=='box' or potential=='well':
        dtau=.025
    else:
        dtau=.01
    print 'dtau should be', dtau
    return dtau

def DefineAnalyticalSoln(potential,center,omega_vibration):
    print 'define analytical soln', omega_vibration
    if potential=='harmonic':
        analyticSoln=(omega_vibration/(np.pi))**(0.25)*np.exp(-(center**2)*omega_vibration/2.0)
        
    elif potential=='box':
        analyticSoln=np.array([ np.sqrt(2./6.)*np.sin(np.pi*(c+3.0)/6.) if cond else 0 for (c,cond) in zip(center,((center>=-3)&(center<=3)))])
    elif potential=='double-well':
        print 'double well!'
        k=omega_vibration**2
        b=omega_vibration
        a=b/2.0
        g=b**2/(4.0*a)
        minLocations=np.sqrt(b/(2.0*a))
        analyticSoln=(omega_vibration/(np.pi))**(0.25)*np.exp(-((center-minLocations)**2)*omega_vibration/2.0)
        analyticSoln=analyticSoln+(omega_vibration/(np.pi))**(0.25)*np.exp(-((center+minLocations)**2)*omega_vibration/2.0)
        print 'yah yah yah'
    elif potential=='well':
        v_n=2.0
        E=.8826
        E=2.0*v_n**2/6.0**2
        k1=np.sqrt(2.0*(2.0-E))
        alph=k1*np.tan(k1*3.0/2.0) 
        analyticSoln=np.array([ np.cos(alph*c) if cond else np.exp(-k1*c) for (c,cond) in zip(center,(center<=3))])
        analyticSoln=np.array([ anaSoln if cond else np.exp(k1*c) for (c,cond,anaSoln) in zip(center,(center>=-3),analyticSoln)])
    elif potential=='halfHO':
        analyticSoln=np.sqrt(2)*(omega_vibration/(np.pi))**(0.25)*center*np.exp(-(center**2)*omega_vibration/2.0)
    else:
        analyticSoln=np.ones(center.size)
    analyticSoln=analyticSoln/np.sqrt(np.sum(analyticSoln*analyticSoln))
    return analyticSoln

def DiffusionMonteCarlo(potentialInfo,N_steps=1000, N=1000, binN=30, alphaVar=.5,omega_vibration=1.0,recrossCorrection=False,PbExp=True):
    print 'Diffusion Omega', omega_vibration
    #Constants and Relationships
    (potential,V,invPotential)=potentialInfo
    dtau=set_dtau(potential,omega_vibration,recrossCorrection)
    alpha=alphaVar/dtau
    reducedMass=1.0
    plotALotOfStuff=False
    plotStuff=False
    plotXV=False

    print 'Initial Parameters \n  potential = ',potential,'\n  alpha = ',alpha,'\n  dtau = ',dtau,'\n  Time steps = ', N_steps, '\n  Number of particles = ',N, '\n  Definition of birth/death',('Exp' if PbExp else 'linear')
    T=N_steps*dtau
    

    D=.25
    sigma_dx=(2*D*dtau)**0.5
    mu_dx=0.0
    
    #non important constants
    colors = np.random.rand(N)
    
    #initial conditions
    mu=0.0
    sigma=.25
    x1,reasonableRange=initialize(potentialInfo, mu, sigma, N,omega_vibration)
    N_size=N
    dx=np.zeros(N_size)

    v=V(x1,omega_vibration)
    v_ref=np.average(v)

    k=v-v_ref  #Which I never actually use

    #Lists of stuff I want to collect
    v_ref_list=[]
    E=[]
    stdevE=[]        
    listOfN_size=[]
    histNormalizedList=[]
    saveStep=0
    hist=np.zeros(binN)
    histNextGen=np.zeros(binN) 
    bin_edges=1
    census=[]
    v_refdataFile=open('dmcVrefData.data','w')
    dataFile=open('dmcData-'+potential+'.data','w')
    dataFile.write('count, step_num    x1             dx              v         v_ref\n')
    count=0
    for step in np.arange(N_steps):
                
        ####### Propagate #######
#        dx=np.random.normal(0,1,N_size)*np.sqrt(dtau)
        #dx is probably a problem for you.  you should look at this again :-( 
        dx=np.random.normal(mu_dx, sigma_dx, N_size)
        
        #if potential=='double-well' or potential=='harmonic': #if it is symmetric... is this cheating
        #swap the sign of the x position of  approximately .1% of the particles
        #   swapers=np.random.random(N_size)\
            #  swapers=np.array([-1 if cond else 1 for cond in (swapers<.1)])
        #  x1=x1*swapers

        x1=x1+dx
        v=V(x1,omega_vibration)#x1,potential,omega=omega_vibration)#-1.0/r1
        k=v-v_ref

        if plotALotOfStuff:
            plt.scatter(x1,dx)
            plt.xlabel('x positions')
            plt.ylabel('dx displacement')
            plt.plot([np.average(x1)]*3,np.arange(-1,2))
            plt.plot(np.arange(-1,2),[np.average(dx)]*3)
            plt.plot([0]*3,np.arange(-1,2),c='k')
            plt.plot(np.arange(-1,2),[0]*3,c='k')
            
            plt.show()
            
            plt.scatter(x1,v)
            print 'vref', v_ref
            plt.show()
            
            print 'max and min',np.max(dx), np.min(dx)

        
        ####### Generate a random number #######
        N_r=np.random.random(N_size)
        
        ######## Deaths #########
        if PbExp:
            P_d=1-np.exp(-(v-v_ref)*dtau)
            #print 'AY YAI YAI'
        else:
            P_d=(v-v_ref)*dtau

        #P_d[P_d < 0] = 0
        Diff=N_r-P_d
        
        #if N_r>P_d, then Diff is positive and then the particle  DOES NOT die and thus survives
        #that's why this mask is called survive
        mask_survive = (Diff>0)  
        tempRecrossCensus=0
        tempCensus=np.sum(np.array(Diff<0).astype(int))
        survivors=x1[mask_survive]

        ###########DEATHS by ReCrossing!################
        if recrossCorrection:
            print 'UHOH'
            P_recrossDeath=np.exp(-2.0*(x1-dx)*x1*np.sqrt(reducedMass*reducedMass)/dtau)
            print 'P_recrossDeath 1', P_recrossDeath
            P_recrossDeath=P_recrossDeath[mask_survive]
            Diff=N_r[mask_survive]-P_recrossDeath  ##Should this be a different random number?
            mask_survive_recross=(Diff>0)
            tempRecrossCensus=np.sum(np.array(Diff<0).astype(int))
            survivors=survivors[mask_survive_recross]            
        else:
            tempRecrossCenus=0

        ####### Regenerate #######
        P_b=-(v-v_ref)*dtau

        #P_b[P_b < 0] = 0
        Diff=N_r-P_b            
        #if Nr is less than Pb,then Diff is less than zero and the particle procreates.
        mask_b = (Diff<0)
        census.append([np.sum(np.array(mask_b).astype(int)), tempCensus])
        next_gen=x1[mask_b]

##OLDVER##        if PbExp:
##OLDVER##            print 'DANGER'
##OLDVER##            P_exp_b=np.exp(-(v-v_ref)*dtau)-1.0
##OLDVER##            weight_P_b=P_exp_b.astype(int)
##OLDVER##            P_b=P_exp_b-weight_P_b
##OLDVER##        else:
##OLDVER##            P_b=-(v-v_ref)*dtau
##OLDVER##            weight_P_b=np.zeros(P_b.size)
##OLDVER##        Diff=N_r-P_b            
##OLDVER##
##OLDVER##        #if Nr is less than Pb,then Diff is less than zero and the particle procreates.
##OLDVER##        mask_b = (Diff<0)
##OLDVER##        #complicated birthing sequences
##OLDVER##        next_gen=x1[mask_b]
##OLDVER##
##OLDVER##        census.append([np.sum(np.array(mask_b).astype(int)), tempCensus, tempRecrossCensus])
##OLDVER##
##OLDVER##        new_pop=next_gen
##OLDVER##        #print new_pop
##OLDVER##        tot=0
##OLDVER##
##OLDVER##        for (particle,weight) in zip(np.array(x1),weight_P_b):
##OLDVER##            if weight>0: #i.e. the dead can't reproduce
##OLDVER##                print 'CAUTION'
##OLDVER##                tot=tot+weight
##OLDVER##                temp=np.tile(particle,weight)
##OLDVER##                new_pop=np.concatenate((new_pop,temp))
##OLDVER##                print 'watch this grow?', new_pop.size
##OLDVER##
##OLDVER##        ######next_gen=new_pop  #maybe not the most efficient memory wise

        #for color coding the birth/death events
        if plotALotOfStuff:# or  step==(N_steps-1):
            colorcoding=['k']*N_size
            for bi,b in enumerate(mask_b):
                if b:
                    colorcoding[bi]='b'
                else:
                    if not mask_survive[bi]:
                        colorcoding[bi]='r'

            plt.xlabel('particle number')
            plt.ylabel('particle x location')
            plt.scatter(np.arange(N_size),x1-dx, c='m',s=N_r*50)
            plt.quiver(np.arange(N_size),x1-dx,0,dx,color=colorcoding)
            plt.plot(np.arange(N_size),[0]*N_size,c='m')
            plt.scatter(np.arange(N_size),x1, c=colorcoding,s=N_r*20)

            ####### Collect Survivors and Next Generation into a set #######
        new_population=np.concatenate((survivors,next_gen))
        N_size=new_population.size
        
        ####### Readjust criteria v_ref #######
        v_average=np.average(V(new_population,omega_vibration))#new_population,potential,omega=omega_vibration))
#        print step, v_average
        v_ref=v_average+(alpha*(1-float(N_size)/float(N)))
        v_ref_list.append(v_ref)

        if plotXV:
            plt.scatter(x1,v)
            plt.pause(.025)
            plt.clf()

            
            
####### Prepare for the next cycle #######
        x1=new_population
#        print 'x1 is at the end of step ', step,' is ', x1
        listOfN_size.append(N_size)    
        #################IF WE ARE HALF WAY THROUGH, SAVE STUFF###################
        if step>=np.floor(N_steps/2):
            v_refdataFile.write(str(step)+'   '+str(v_ref)+'\n')
            for guy in zip(x1,V(x1,omega_vibration)):
                count=count+1
                dataFile.write(str(str(count)+'       '+str(step)+'       '+str(guy[0])+' '+str(guy[1])+' '+str(v_ref)+'\n'))
            dataFile.write('\n')
            histTempNextGen,bin_edges=np.histogram(next_gen,bins=binN,range=reasonableRange)
            histTemp,bin_edges=np.histogram(x1,bins=binN,range=reasonableRange)
            hist=hist+histTemp
            histNextGen=histTempNextGen+histNextGen
        E.append(v_ref)
        stdevE.append(np.std(v))
        if step==np.floor(N_steps/2):
            print 'half way data!', np.max(x1), np.min(x1), np.average(x1), np.std(x1)
        saveStep=step
        if N_size<1:
            print 'MASSIVE DIE OFF'
            print step
            break
        
    print 'census data', step, census[step],'vref',v_ref, 'Nsize',N_size

    iterations=saveStep
    #normalize histogram at the end of the day
        
    width=0.7*(bin_edges[1]-bin_edges[0])
    center=(bin_edges[:-1]+bin_edges[1:])/2.0
    histPlot=hist/np.sqrt(np.sum(hist*hist))

    analyticSoln=DefineAnalyticalSoln(potential,center,omega_vibration) #(omega_vibration/(1))**(0.25)*np.exp(-(center**2)*omega_vibration/2)
    expectationHist=histNextGen
    
    hist=hist/np.sqrt(np.sum(hist*hist))
    
    print '\n Final Results: \n   Average E = ',np.average(E), '\n   Standard Deviation of E = ',np.std(E), '\n   Final Population =' ,N_size ,'\n   Cenus (birth, death, recrossDeath?) = ', np.average(np.array(census),axis=0), '+/-',np.std(np.array(census),axis=0),'\n\n'
    return iterations, E, N_size, hist,expectationHist,bin_edges, np.array(histNormalizedList)

def multipledmc(multiple,potentialInfo,omega=1.0,recrossing=False,ExpBirth=True):

	totalAveE=[]
        (potential,Vpotential,invPotential)=potentialInfo
        print 'omega multiple dmc', omega
	totalN=[]
	totalIterations=[]
	totalStdE=[]
	nboxes=101
	hist=np.zeros(nboxes)
        nextGenHist=np.zeros(nboxes)
	timesteps=2000
        N_initial=100
        a=.5**2
	for j in range(0,multiple):
	    print 'Simulation Number ', j
	    iterations,E,N,histTemp,nextHistTemp,bin_edges,histNormalizedList=DiffusionMonteCarlo(potentialInfo,N_steps=timesteps,N=N_initial,binN=nboxes,alphaVar=a,omega_vibration=omega,recrossCorrection=recrossing,PbExp=ExpBirth)

            plt.plot(E)
            plt.show()
	    totalAveE.append(np.average(E))

            hist=hist+histTemp
	    nextGenHist=nextGenHist+nextHistTemp
	    #histTemp=hist/np.sqrt(np.sum(hist*hist))
	    width=0.7*(bin_edges[1]-bin_edges[0])
            center=(bin_edges[:-1]+bin_edges[1:])/2.0
            np.savetxt('run'+str(j)+'amplitude.data',zip(center,histTemp))
            print 'birth Hist Temp'
	    plt.bar(center,nextHistTemp,align='center',width=width)
            plt.show()

            print 'histTemp'
            plt.bar(center,histTemp/np.sqrt(np.sum(histTemp*histTemp)),align='center',width=width)

	    analyticSoln=DefineAnalyticalSoln(potential,center,omega) 
	    plt.plot(center,analyticSoln,'.r-')
	    plt.savefig('HistogramExpectation'+str(j)+'.pdf')
            plt.show()
#	    plt.pause(.5)
#	    plt.clf()
            
	    
	    totalStdE.append(np.std(E))
	    totalN.append(N)
	    totalIterations.append(iterations)
	    #print 'Simulation',j,'  Total iterations',iterations, '  Size at end', N
	hist=hist/np.sqrt(np.sum(hist*hist))
        nextGenHist=nextGenHist/np.sqrt(np.sum(nextGenHist*nextGenHist))
        
        width=0.7*(bin_edges[1]-bin_edges[0])
        center=(bin_edges[:-1]+bin_edges[1:])/2.0
        print 'histTemp'
        plt.bar(center,hist,align='center',width=width)
        plt.show()
        print 'birth Hist Temp'
        plt.bar(center,nextGenHist,align='center',width=width)
        plt.show()

        expectationX=np.sum(hist*nextGenHist*center)
        
        secondMoment=np.sum(hist*nextGenHist*center*center)
        variance=secondMoment-expectationX**2
        print 'expectation X' , expectationX, 'second moment', secondMoment, 'variance', variance
	print '\n Average Energy From all Simulations: ', np.average(np.array(totalAveE)),'+/-',np.std(np.array(totalAveE))    
	print '\n Average number of steps:', np.average(np.array(totalIterations)),'+/-',np.std(np.array(totalIterations))
	#fig=plt.figure()
	#ax1=fig.add_axes([0.1,0.55,0.8,0.35])
	#ax2=fig.add_axes([0.1,0.1,0.8,0.35])
	#ax2.set_xlabel('x position')
	#
	##ax3=ax2.twinx()
	##ax1.set_ylabel('Total Iterations Completed',color='r')
	#ax1.set_ylabel('Number of particles',color='k')
	#ax1.set_xlabel('Simulation Number')
	##ax1.set_xtics(np.arange(1000,1000+j*500,500))
	##ax1.plot(np.arange(len(totalIterations)),np.array(totalIterations),'.r-')
	#ax1.plot(np.arange(len(totalN)),np.array(totalN), '.k-')
	#ax1b=ax1.twinx()
	#ax1b.set_ylabel('Average Energy (E_ref) from the final N-100 steps',color='b')
	#ax1b.errorbar(np.arange(len(totalAveE)),np.array(totalAveE),np.std(totalAveE),c='b')
	#
	#
	histNorm=hist/np.sqrt(np.sum(hist*hist))

	width=0.7*(bin_edges[1]-bin_edges[0])
	center=(bin_edges[:-1]+bin_edges[1:])/2.0
	#ax2.bar(center,histNorm,align='center',width=width)
	analyticSoln=DefineAnalyticalSoln(potential,center,omega)
        #
        #print 'center ',center[0:3],center[-3:]
	#
	#ax2.plot(center,analyticSoln,'.r-')
	#plt.savefig('Energy_and_Histogram.pdf')
	#plt.show()

        return histNorm, center, np.average(np.array(totalAveE)), np.std(np.array(totalAveE))

print 'done with dmc'
