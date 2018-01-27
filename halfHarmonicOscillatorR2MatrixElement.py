import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import os
au2wn=219474.63
nWalkers=1000


Wfn1Left =dmc.wavefunction(nWalkers,'half harmonic left',plotting=False)
Wfn1Right=dmc.wavefunction(nWalkers,'half harmonic right',plotting=False)
idealVref=Wfn1Left.getIdealVref()

nReps=5
nRepsDesc=5
nSteps=100
nDesSteps=25
v_ref_Left=[]
pop_Left=[]
v_ref_array_Left=np.zeros((nReps,nSteps))
v_ref_desc_array_Left=np.zeros((2,nReps,nDesSteps))
pop_array_Left=np.zeros((nReps,nSteps))
pop_desc_array_Left=np.zeros((2,nReps,nDesSteps))
v_ref_Right=[]
pop_Right=[]
v_ref_array_Right=np.zeros((nReps,nSteps))
v_ref_desc_array_Right=np.zeros((2,nReps,nDesSteps))
pop_array_Right=np.zeros((nReps,nSteps))
pop_desc_array_Right=np.zeros((2,nReps,nDesSteps))
nBins=51

Destination='State1-nT'+str(nDesSteps)+'nR'+str(nReps)+'nRD'+str(nRepsDesc)+'Pop'+str(nWalkers)+'/'
print 'Destination:', Destination
file_path=Destination+'logFile.data'
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
        os.makedirs(directory)

xhists_1=np.zeros((3,nReps,nBins)) #Psi, Psi_popPsi_pop, Psi_popPsi_vref
AveDescendantWeightsLeftPop=[]
AveDescendantWeightsLeftVRef=[]
AveDescendantWeightsRightPop=[]
AveDescendantWeightsRightVRef=[]
R2_LeftPop=np.zeros((nReps))
R2_LeftVRef=np.zeros((nReps))
R2_RightPop=np.zeros((nReps))
R2_RightVRef=np.zeros((nReps))

Wfn1Left.setX(Wfn1Left.xcoords-.5)
Wfn1Right.setX(Wfn1Right.xcoords+.5)
for n in range(nReps):
    vref,pop,x,d=Wfn1Left.propagate(Wfn1Left.xcoords,nSteps)
    vref=np.array(vref)
    v_ref_Left.append(np.average(vref[nSteps/2:])*au2wn)
    pop_Left.append(np.average(pop[nSteps/2:]))
    v_ref_array_Left[n,:]=vref*au2wn
    pop_array_Left[n,:]=pop
    Psi1LeftHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True)
    desc_weight_array=np.zeros((nRepsDesc,x.size))
    #for psi squared
    #first the usual way, adjusting vref during simulation
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Left.propagate(x,nDesSteps,nSize=nWalkers)
        v_ref_desc_array_Left[0,n,:]=np.array(vrd)*au2wn
        pop_desc_array_Left[0,n,:]=popd
        desc_weight_array[i]=descendantWeights
    AveDescendantWeightsLeftPop.append(np.average(desc_weight_array,axis=0))
    Psi2LeftHist_pop_pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=AveDescendantWeightsLeftPop[n])
    #second the unusual way, setting vref to be the 'proper' value
    desc_weight_array=np.zeros((nRepsDesc,x.size))
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Left.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref)
        v_ref_desc_array_Left[1,n,:]=np.array(vrd)*au2wn
        pop_desc_array_Left[1,n,:]=popd
        desc_weight_array[i]=descendantWeights

    AveDescendantWeightsLeftVRef.append(np.average(desc_weight_array,axis=0))
    Psi2LeftHist_pop_vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=AveDescendantWeightsLeftVRef[n])
    R2_LeftPop[n]=np.sum(x*x*AveDescendantWeightsLeftPop[n])/np.sum(AveDescendantWeightsLeftPop[n])
    R2_LeftVRef[n]=np.sum(x*x*AveDescendantWeightsLeftVRef[n])/np.sum(AveDescendantWeightsLeftVRef[n])
##########
##########
    vref,pop,x,d=Wfn1Right.propagate(Wfn1Right.xcoords,nSteps)
    vref=np.array(vref)
    v_ref_Right.append(np.average(vref[nSteps/2:])*au2wn)
    pop_Right.append(np.average(pop[nSteps/2:]))
    v_ref_array_Right[n,:]=vref*au2wn
    pop_array_Right[n,:]=pop
    Psi1RightHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True)
    desc_weight_array=np.zeros((nRepsDesc,x.size))
    #for psi squared
    #first the usual way, adjusting vref during simulation
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Right.propagate(x,nDesSteps,nSize=nWalkers)
        v_ref_desc_array_Right[0,n,:]=np.array(vrd)*au2wn
        pop_desc_array_Right[0,n,:]=popd
        desc_weight_array[i]=descendantWeights
    AveDescendantWeightsRightPop.append(np.average(desc_weight_array,axis=0))
    Psi2RightHist_pop_pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=AveDescendantWeightsRightPop[n])
    #second the unusual way, setting vref to be the 'proper' value
    desc_weight_array=np.zeros((nRepsDesc,x.size))
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Right.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref)
        v_ref_desc_array_Right[1,n,:]=np.array(vrd)*au2wn
        pop_desc_array_Right[1,n,:]=popd
        desc_weight_array[i]=descendantWeights

    AveDescendantWeightsRightVRef.append(np.average(desc_weight_array,axis=0))
    Psi2RightHist_pop_vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=AveDescendantWeightsRightVRef[n])
    
    R2_RightPop[n]=np.sum(x*x*AveDescendantWeightsRightPop[n])/np.sum(AveDescendantWeightsRightPop[n])
    R2_RightVRef[n]=np.sum(x*x*AveDescendantWeightsRightVRef[n])/np.sum(AveDescendantWeightsRightVRef[n])


    xhists_1[0,n]=Psi1LeftHist-Psi1RightHist
    xhists_1[1,n]=Psi2LeftHist_pop_pop+Psi2RightHist_pop_pop
    xhists_1[2,n]=Psi2LeftHist_pop_vref+Psi2RightHist_pop_vref
    

print 'average energy',np.average(v_ref_Left)
print 'standard deviation', np.std(v_ref_Left)
print 'uncertainity', (np.max(v_ref_Left)-np.min(v_ref_Left))/(2.0*np.sqrt(nReps))
print 'average population', np.average(pop_Left)

print 'Ideal is:', idealVref*au2wn

## Calculate V_1|R2|V1
print 'Analytical R^2 is', 1.5*1.0/Wfn1Left.getAlpha()
print 'Numerical', R2_LeftPop,R2_LeftVRef, R2_RightPop,R2_RightVRef
print 'averaged', np.average(R2_LeftPop+R2_RightPop), np.average(R2_LeftVRef+R2_RightVRef)
colors=['r','b','k','g']
#plotting first v_Ref
for n in range(nReps):
    plt.plot(np.arange(nSteps),v_ref_array_Left[n],c='r')
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_Left[0,n],c='m')
plt.show()    
for n  in range(nReps):
    plt.plot(np.arange(nSteps),pop_array_Left[n],c='r')
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_Left[0,n],c='m')
plt.show()

bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
for n in range(nReps):
    plt.plot(bin_center,xhists_1[0,n],c='r')
plt.show()

for n in range(nReps):
    plt.plot(bin_center,xhists_1[1,n],c='r')
    plt.plot(bin_center,xhists_1[2,n],c='m')
plt.show()

print 'done!'

