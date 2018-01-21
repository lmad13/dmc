import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
nWalkers=2000
Wfn0=dmc.wavefunction(nWalkers,'harmonic',plotting=False)
v_ref=[]
pop_0=[]
nReps=25
nSteps=1000

colors=['r','b','k','g']
v_ref_array_0=np.zeros((nReps,nSteps))
pop_array_0=np.zeros((nReps,nSteps))
xhists_0=np.zeros((nReps,50))

for n in range(nReps):
    vref,pop,x,d=Wfn0.propagate(Wfn0.xcoords,nSteps)
    vref=np.array(vref)
    plt.plot(np.arange(nSteps),vref*au2wn,c='r')
    v_ref.append(np.average(vref[nSteps/2:])*au2wn)
    pop_0.append(np.average(pop[nSteps/2:]))
    v_ref_array[n,:]=vref*au2wn
    pop_array[n,:]=pop
print 'average energy',np.average(v_ref)
print 'standard deviation', np.std(v_ref)
print 'uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print 'average population', np.average(pop_0)
np.savetxt('vrefForNOTFixed.data', v_ref_array.T)
np.savetxt('popForNOTFixed.data', pop_array.T)
print Wfn0.omega/2.0
plt.show()
nReps=25
Wfn1=dmc.wavefunction(nWalkers,'harmonic',plotting=False)
v_ref_1=[]
pop_1=[]
v_ref_array=np.zeros((nReps,nSteps))
pop_array=np.zeros((nReps,nSteps))
for n in range(nReps):
    vref,pop,x,d=Wfn1.propagate(Wfn1.xcoords,nSteps,setV_ref=True,ConstantV_ref=1000.0/au2wn)
    v_ref_1.append(np.average(vref[nSteps/2:])*au2wn)
    v_ref_array[n,:]=vref
    pop_array[n,:]=pop
    plt.plot(np.arange(nSteps),pop,c='b')
    pop_1.append(np.average(pop[nSteps/2:]))
print 'v1 average energy',np.average(v_ref_1)
print 'v1 standard deviation', np.std(v_ref_1)
print 'v1 uncertainity', (np.max(v_ref_1)-np.min(v_ref_1))/(2.0*np.sqrt(nReps))
print 'pop', np.average(pop_1)
print 'v1 standard deviation', np.std(pop_1)
print 'v1 uncertainity', (np.max(pop_1)-np.min(pop_1))/(2.0*np.sqrt(nReps))
np.savetxt('popForVRefFIXED.data', pop_array.T)

print 1.0*Wfn1.omega/2.0
plt.show()

#wavefunction class

#initalize wavefunction
#place walkers on HO

#function propagate
#propagate wavefunction for nsteps WITHOUT kmnowing the final v_ref
#take in walker positions, V, nsteps, delta t
#return vref measurment (avearged over 1000 steps, only at equilibrium), snapshots of wavefunction every 500 t steps after equilibration

#function CalcDesc
#calculate descendant weights
#take in walker positions, V, nsteps, delta t         
#return the number of descendants



