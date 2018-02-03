#!/usr/bin/python

import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
au2wn=219474.63

if len(sys.argv)<5:
	print 'usage: ./scipt.py nWalkers nReps nSteps nDesSteps nRepsDesc'
	end
else:
	print sys.argv

nWalkers=int(sys.argv[1])
nReps=int(sys.argv[2])
nSteps=int(sys.argv[3])
nDesSteps=int(sys.argv[4])
nRepsDesc=int(sys.argv[5])

Wfn0 =dmc.wavefunction(nWalkers,'harmonic',plotting=False)

idealVref=Wfn0.getIdealVref()

v_ref=[]
pop0=[]
v_ref_array=np.zeros((nReps,nSteps))
v_ref_desc_array=np.zeros((2,nReps,nDesSteps))
pop_array=np.zeros((nReps,nSteps))
pop_desc_array=np.zeros((2,nReps,nDesSteps))
nBins=51

Destination='State0-nT'+str(nDesSteps)+'nR'+str(nReps)+'nRD'+str(nRepsDesc)+'Pop'+str(nWalkers)+'/'
print 'Destination:', Destination
file_path=Destination+'logFile.data'
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
        os.makedirs(directory)

xhists_0=np.zeros((3,nReps,nBins)) #Psi, Psi_popPsi_pop, Psi_popPsi_vref
AveDescendantWeightsPop=[]
AveDescendantWeightsVRef=[]

R2_Pop=np.zeros((nReps))
R2_VRef=np.zeros((nReps))


for n in range(nReps):
    vref,pop,x,d=Wfn0.propagate(Wfn0.xcoords,nSteps)

    Psi0Hist,bin_edges=np.histogram(x, bins=nBins, range=(-2.0,2.0),density=True)
    vref=np.array(vref)
    v_ref.append(np.average(vref[nSteps/2:])*au2wn)
    pop0.append(np.average(pop[nSteps/2:]))
    v_ref_array[n,:]=vref*au2wn
    pop_array[n,:]=pop

    desc_weight_array=np.zeros((nRepsDesc,x.size))
    #for psi squared
    #first the usual way, adjusting vref during simulation
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn0.propagate(x,nDesSteps,nSize=nWalkers)
        v_ref_desc_array[0,n,:]=np.array(vrd)*au2wn
        pop_desc_array[0,n,:]=popd
        desc_weight_array[i]=descendantWeights
    AveDescendantWeightsPop.append(np.average(desc_weight_array,axis=0))
    Psi2Hist_pop_pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=AveDescendantWeightsPop[n])
    #second the unusual way, setting vref to be the 'proper' value
    desc_weight_array=np.zeros((nRepsDesc,x.size))
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn0.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref)
        v_ref_desc_array[1,n,:]=np.array(vrd)*au2wn
        pop_desc_array[1,n,:]=popd
        desc_weight_array[i]=descendantWeights

    AveDescendantWeightsVRef.append(np.average(desc_weight_array,axis=0))
    Psi2Hist_pop_vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True,weights=AveDescendantWeightsVRef[n])
    R2_Pop[n]=np.sum(x*x*AveDescendantWeightsPop[n])/np.sum(AveDescendantWeightsPop[n])
    R2_VRef[n]=np.sum(x*x*AveDescendantWeightsVRef[n])/np.sum(AveDescendantWeightsVRef[n])

##########

    xhists_0[0,n]=Psi0Hist
    xhists_0[1,n]=Psi2Hist_pop_pop
    xhists_0[2,n]=Psi2Hist_pop_vref
    

print 'average energy',np.average(v_ref)
print 'standard deviation', np.std(v_ref)
print 'uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print 'average population', np.average(pop0)

print 'Ideal is:', idealVref*au2wn

## Calculate V_1|R2|V1
print 'Analytical R^2 is', 0.5/Wfn0.getAlpha()
#print 'Numerical', R2_LeftPop,R2_LeftVRef, R2_RightPop,R2_RightVRef
print 'averaged', np.average(R2_Pop), np.average(R2_VRef)
print 'standard deviation', np.std(R2_Pop),np.std(R2_VRef)
unc=[(np.max(R2_Pop)-np.min(R2_Pop))/(2.0*np.sqrt(nReps)),(np.max(R2_VRef)-np.min(R2_VRef))/(2.0*np.sqrt(nReps))]
print 'uncertainity', unc
fileR2Data=open('R2Data-GroundState.data','a')
fileR2Data.write(str(0)+'   '+str(nWalkers)+'   '+str(nSteps)+'   '+str(nReps)+'   '+str(nDesSteps)+'   '+str(nRepsDesc)+'      ')
fileR2Data.write(str(0.5*1.0/Wfn0.getAlpha())+'   '+str(np.average(R2_Pop))+'   '+str(np.average(R2_VRef))+'   '+str(np.std(R2_Pop))+'   '+str(np.std(R2_VRef))+'   '+str(unc[0])+'   '+str(unc[1])+'\n')
fileR2Data.close()
colors=['r','b','k','g']
red='#f60000'
darkRed='#9c0202'
orange='#ff9a00'
darkOrange='#d86200'
purple='#ad4dff'
darkPurple='#5a00b8'
pink='#fe7eff'
darkPink='#b900ba'
blue='#4d62ff'
darkBlue='#0c008d'
green='#01b24c'
darkGreen='#005423'
#plotting first v_Ref
for n in range(nReps):
    plt.plot(np.arange(nSteps),v_ref_array[n],c='blue')
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array[0,n],c='pink')
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array[1,n],c='darkgreen')
plt.ylabel('Energy (1/cm)')
plt.xlabel('Step ('+str(Wfn0.get_dtau())+')')
plt.savefig(Destination+'Step-Energy.png')
plt.clf()
for n  in range(nReps):
    plt.plot(np.arange(nSteps),pop_array[n],c='blue')
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array[0,n],c='pink')
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array[1,n],c='darkgreen')
plt.ylabel('Population (Walkers)')
plt.xlabel('Step ('+str(Wfn0.get_dtau())+')')
plt.savefig(Destination+'Step-Pop.png')
plt.clf()

bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
for n in range(nReps):
    plt.plot(bin_center,xhists_0[0,n],c='r')
plt.ylabel('Amplitude (Psi)')
plt.xlabel('X (bohr)')
plt.savefig(Destination+'Psi_0.png')
plt.clf()

for n in range(nReps):
    plt.plot(bin_center,xhists_0[1,n],c='r')
    plt.plot(bin_center,xhists_0[2,n],c='m')
plt.ylabel('Amplitude (Psi*Psi)')
plt.xlabel('X (bohr)')
plt.savefig(Destination+'Psi_0Psi_0.png')
plt.clf()

print 'done!'

