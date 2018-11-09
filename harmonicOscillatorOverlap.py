import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
au2wn=219474.63
nWalkers=4000


Wfn0=dmc.wavefunction(nWalkers,'harmonic',plotting=False)
idealVref=Wfn0.getIdealVref()
v_ref=[]
pop_0=[]
nReps=25
nSteps=1000
nDesSteps=25
v_ref_array_0=np.zeros((nReps,nSteps))
v_ref_desc_array_0=np.zeros((2,nReps,nDesSteps))
pop_array_0=np.zeros((nReps,nSteps))
pop_desc_array_0=np.zeros((2,nReps,nDesSteps))
nBins=101

xhists_0=np.zeros((3,nReps,nBins)) #Psi, Psi_popPsi_pop, Psi_popPsi_vref

Wfn0.setX(Wfn0.xcoords+.1)
for n in range(nReps):
    vref,pop,x,d=Wfn0.propagate(Wfn0.xcoords,nSteps)
    vref=np.array(vref)
    v_ref.append(np.average(vref[nSteps/2:])*au2wn)
    pop_0.append(np.average(pop[nSteps/2:]))
    v_ref_array_0[n,:]=vref*au2wn
    pop_array_0[n,:]=pop
    PsiHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True)

    #for psi squared
    #first the usual way, adjusting vref during simulation
    vrd,popd,xd,descendantWeights=Wfn0.propagate(x,nDesSteps,nSize=nWalkers)
    v_ref_desc_array_0[0,n,:]=np.array(vrd)*au2wn
    pop_desc_array_0[0,n,:]=popd
    Psi2Hist_pop_pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=descendantWeights)
    #second the unusual way, setting vref to be the 'proper' value
    vrd,popd,xd,descendantWeights=Wfn0.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref)
    ### check herer for UW plot agreement
    v_ref_desc_array_0[1,n,:]=np.array(vrd)*au2wn
    pop_desc_array_0[1,n,:]=popd
    Psi2Hist_pop_vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=descendantWeights)

    xhists_0[0,n]=PsiHist
    xhists_0[1,n]=Psi2Hist_pop_pop
    xhists_0[2,n]=Psi2Hist_pop_vref
    

print 'average energy',np.average(v_ref)
print 'standard deviation', np.std(v_ref)
print 'uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print 'average population', np.average(pop_0)
np.savetxt('GroundStateData/GroundStateData/vrefForPOP.data', v_ref_array_0.T)
np.savetxt('GroundStateData/popForPOP.data', pop_array_0.T)
np.savetxt('GroundStateData/vrefDescFor-POP-POP.data',v_ref_desc_array_0[0].T)
np.savetxt('GroundStateData/vrefDescFor-POP-VREF.data',v_ref_desc_array_0[1].T)
np.savetxt('GroundStateData/popDescFor-POP-VREF.data',pop_desc_array_0[1].T)
np.savetxt('GroundStateData/popDescFor-POP-POP.data',pop_desc_array_0[0].T)

print 'Ideal is:', idealVref*au2wn

Wfn1=dmc.wavefunction(nWalkers,'harmonic',plotting=False)
Wfn1.setX(Wfn1.xcoords+.1)
idealVref=Wfn0.getIdealVref()
v_ref_1=[]
pop_1=[]
v_ref_array_1=np.zeros((nReps,nSteps))
v_ref_desc_array_1=np.zeros((2,nReps,nDesSteps))
pop_desc_array_1=np.zeros((2,nReps,nDesSteps))
pop_array_1=np.zeros((nReps,nSteps))
xhists_1=np.zeros((3,nReps,nBins))
for n in range(nReps):
    vref,pop,x,d=Wfn1.propagate(Wfn1.xcoords,nSteps,setV_ref=True,ConstantV_ref=idealVref)
    vref=np.array(vref)
    v_ref_1.append(np.average(vref[nSteps/2:])*au2wn)
    v_ref_array_1[n,:]=vref*au2wn
    pop_array_1[n,:]=pop
    pop_1.append(np.average(pop[nSteps/2:]))
    PsiHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True)
    xhists_1[0,n]=PsiHist
    
    vrd,popd,xd,descendantWeights=Wfn1.propagate(x,nDesSteps,setV_ref=False,nSize=x.size)
    v_ref_desc_array_1[0,n,:]=np.array(vrd)*au2wn
    pop_desc_array_1[0,n,:]=popd
    Psi2_vref_popHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=descendantWeights)
    xhists_1[1,n]=Psi2_vref_popHist

    vrd,popd,xd,descendantWeights=Wfn1.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref)
    v_ref_desc_array_1[1,n,:]=np.array(vrd)*au2wn
    pop_desc_array_1[1,n,:]=popd
    Psi2_vref_vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=descendantWeights)
    xhists_1[2,n]=Psi2_vref_vref

print 'v1 average energy',np.average(v_ref_1)
print 'v1 standard deviation', np.std(v_ref_1)
print 'v1 uncertainity', (np.max(v_ref_1)-np.min(v_ref_1))/(2.0*np.sqrt(nReps))
print 'pop', np.average(pop_1)
print 'v1 standard deviation', np.std(pop_1)
print 'v1 uncertainity', (np.max(pop_1)-np.min(pop_1))/(2.0*np.sqrt(nReps))
np.savetxt('GroundStateData/popForVRefFIXED.data', pop_array_1.T)
np.savetxt('GroundStateData/vref-VREF-descendants-POP.data',v_ref_desc_array_1[0].T)
np.savetxt('GroundStateData/popDescFor-POP-POP.data',pop_desc_array_1[0].T)
np.savetxt('GroundStateData/popDescFor-POP-VREF.data.data',pop_desc_array_1[1].T)

print 'ideal is',idealVref*au2wn

colors=['r','b','k','g']
#plotting first v_Ref
for n in range(nReps):
    plt.plot(np.arange(nSteps),v_ref_array_0[n],c='r')
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[0,n],c='m')
    plt.plot(np.arange(nSteps),v_ref_array_1[n],c='b')
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_1[0,n],c='g')
plt.show()    
for n  in range(nReps):
    plt.plot(np.arange(nSteps),pop_array_0[n],c='r')
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[0,n],c='m')
    plt.plot(np.arange(nSteps),pop_array_1[n],c='b')
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_1[0,n],c='g')
plt.show()

bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
for n in range(nReps):
    plt.plot(bin_center,xhists_0[0,n],c='r')
    plt.plot(bin_center,xhists_1[0,n],c='b')
plt.show()

for n in range(nReps):
    plt.plot(bin_center,xhists_0[1,n],c='r')
    plt.plot(bin_center,xhists_0[2,n],c='m')
    plt.plot(bin_center,xhists_1[1,n],c='b')
    plt.plot(bin_center,xhists_1[2,n],c='g')
plt.show()

np.savetxt('GroundStateData/PsiHist-pop.data',xhists_0[0].T)
np.savetxt('GroundStateData/PsiHist-AVERAGE-pop.data',zip(bin_center,np.average(xhists_0[0],axis=0)))
np.savetxt('GroundStateData/Psi2Hist-pop-pop.data',xhists_0[1].T)
np.savetxt('GroundStateData/Psi2Hist-AVERAGE-pop-pop.data',zip(bin_center,np.average(xhists_0[1],axis=0)))
np.savetxt('GroundStateData/Psi2Hist-pop-vref.data',xhists_0[2].T)
np.savetxt('GroundStateData/Psi2Hist-AVERAGE-pop-vref.data',zip(bin_center,np.average(xhists_0[2],axis=0)))

np.savetxt('GroundStateData/PsiHist-vref.data',xhists_1[0].T)
np.savetxt('GroundStateData/PsiHist-AVERAGE-vref.data',zip(bin_center,np.average(xhists_1[0],axis=0)))
np.savetxt('GroundStateData/Psi2Hist-vref-pop.data',xhists_1[1].T)
np.savetxt('GroundStateData/Psi2Hist-AVERAGE-vref-pop.data',zip(bin_center,np.average(xhists_1[1],axis=0)))
np.savetxt('GroundStateData/Psi2Hist-vref-vref.data',xhists_1[2].T)
np.savetxt('GroundStateData/Psi2Hist-AVERAGE-vref-vref.data',zip(bin_center,np.average(xhists_1[2],axis=0)))

print 'done!'

