import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
nWalkers=2000

Wfn0=dmc.wavefunction(nWalkers,'half harmonic',plotting=False)
v_ref=[]
pop_0=[]
nReps=5
nSteps=1000

v_ref_array_0=np.zeros((nReps,nSteps))
pop_array_0=np.zeros((nReps,nSteps))
nBins=30
xhists_0=np.zeros((nReps,nBins))

Wfn0.setX(Wfn0.xcoords+.001)
for n in range(nReps):
    vref,pop,x,d=Wfn0.propagate(Wfn0.xcoords,nSteps)
    vref=np.array(vref)
    v_ref.append(np.average(vref[nSteps/2:])*au2wn)
    pop_0.append(np.average(pop[nSteps/2:]))
    v_ref_array_0[n,:]=vref*au2wn
    pop_array_0[n,:]=pop
    vrd,popd,xd,descendantWeights=Wfn0.propagate(x,25)
    xHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=descendantWeights)
    xhists_0[n]=xHist

print 'average energy',np.average(v_ref)
print 'standard deviation', np.std(v_ref)
print 'uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print 'average population', np.average(pop_0)
np.savetxt('vrefForNOTFixed.data', v_ref_array_0.T)
np.savetxt('popForNOTFixed.data', pop_array_0.T)
print Wfn0.omega/2.0

Wfn1=dmc.wavefunction(nWalkers,'half harmonic',plotting=False)
Wfn1.setX(Wfn1.xcoords+.1)
v_ref_1=[]
pop_1=[]
v_ref_array_1=np.zeros((nReps,nSteps))
pop_array_1=np.zeros((nReps,nSteps))
xhists_1=np.zeros((nReps,nBins))
for n in range(nReps):
    vref,pop,x,d=Wfn1.propagate(Wfn1.xcoords,nSteps,setV_ref=True,ConstantV_ref=3000.0/au2wn)
    vref=np.array(vref)
    v_ref_1.append(np.average(vref[nSteps/2:])*au2wn)
    v_ref_array_1[n,:]=vref*au2wn
    pop_array_1[n,:]=pop
    pop_1.append(np.average(pop[nSteps/2:]))
    vrd,popd,xd,descendantWeights=Wfn1.propagate(x,25,setV_ref=True,ConstantV_ref=3000.0/au2wn)
    xHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=descendantWeights)
    xhists_1[n]=xHist
print 'v1 average energy',np.average(v_ref_1)
print 'v1 standard deviation', np.std(v_ref_1)
print 'v1 uncertainity', (np.max(v_ref_1)-np.min(v_ref_1))/(2.0*np.sqrt(nReps))
print 'pop', np.average(pop_1)
print 'v1 standard deviation', np.std(pop_1)
print 'v1 uncertainity', (np.max(pop_1)-np.min(pop_1))/(2.0*np.sqrt(nReps))
np.savetxt('popForVRefFIXED.data', pop_array_1.T)

print 1.0*Wfn1.omega/2.0

colors=['r','b','k','g']
#plotting first v_Ref
for n in range(nReps):
    plt.plot(np.arange(nSteps),v_ref_array_0[n],c='r')
    plt.plot(np.arange(nSteps),v_ref_array_1[n],c='b')
plt.show()    
for n  in range(nReps):
    plt.plot(np.arange(nSteps),pop_array_0[n],c='r')
    plt.plot(np.arange(nSteps),pop_array_1[n],c='b')
plt.show()

bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
for n in range(nReps):
    plt.plot(bin_center,xhists_0[n],c='r')
    plt.plot(bin_center,xhists_1[n],c='b')
plt.show()



