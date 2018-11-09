#!/usr/bin/python

import os
import sys
import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
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


Wfn0=dmc.wavefunction(nWalkers,'harmonic',plotting=False)
idealVref_0=Wfn0.getIdealVref()

Wfn1Left=dmc.wavefunction(nWalkers/2,'half harmonic left',plotting=False)
Wfn1Left.setX(Wfn1Left.xcoords-.1)
idealVref_1=Wfn1Left.getIdealVref()

Wfn1Right=dmc.wavefunction(nWalkers/2,'half harmonic right',plotting=False)
Wfn1Right.setX(Wfn1Right.xcoords+.1)
idealVref_1=Wfn1Right.getIdealVref()

v_ref=[]
pop_0=[]


v_ref_array_0=np.zeros((nReps,nSteps))
v_ref_desc_array_0=np.zeros((3,2,nReps,nDesSteps))
pop_array_0=np.zeros((nReps,nSteps))
pop_desc_array_0=np.zeros((3,2,nReps,nDesSteps))
xCoord0=[]
R2_0Pop=np.zeros((nReps))   #<0|R2|0>
R2_0VRef=np.zeros((nReps))  #<0|R2|0>
R2_1Pop=np.zeros((nReps))   #<1|R2|1>
R2_1VRef=np.zeros((nReps))  #<1|R2|1>
nBins=51
AveDescendantWeightsGroundPop=[] #will be indexed nReps, x.size  d^0_pop
AveDescendantWeightsGroundVRef=[] #will be indexed nReps, x.size d^0_vref
AveDescendantWeightsExcLeftPop=[] #will be indexed nReps, x.size 
AveDescendantWeightsExcLeftVRef=[] #will be indexed nReps, x.size
AveDescendantWeightsExcRightPop=[] #will be indexed nReps, x.size
AveDescendantWeightsExcRightVRef=[] #will be indexed nReps, x.size
AveDescendantWeightsExcPop=[] #will be indexed nReps, x.size  d^1_pop
AveDescendantWeightsExcVRef=[] #will be indexed nReps, x.size d^1_vref

xhists_0=np.zeros((5,nReps,nBins)) #Psi, Psi_popPsi_pop, Psi_popPsi_vref

Wfn0.setX(Wfn0.xcoords+.1)
##calculate expectation values
Destination='States0to1-nT'+str(nDesSteps)+'nR'+str(nReps)+'nRD'+str(nRepsDesc)+'Pop'+str(nWalkers)+'/'
print 'Destination:', Destination
file_path=Destination+'logFile.data'
directory = os.path.dirname(file_path)
if not os.path.exists(directory):
        os.makedirs(directory)
fileLog=open(file_path,'w')
fileLog.write('nReps       '+str(nReps)+'\n')
fileLog.write('nRepsDesc   '+str(nRepsDesc)+'\n')
fileLog.write('nSteps      '+str(nSteps)+'\n')
fileLog.write('nDesSteps   '+str(nDesSteps)+'\n')
fileLog.write('pop         '+str(nWalkers)+'\n')
fileLog.close()

for n in range(nReps):
    vref,pop,x,d=Wfn0.propagate(Wfn0.xcoords,nSteps)
    xCoord0.append(x)
    desc_weight_array_0=np.zeros((nRepsDesc,x.size))  
    AveDescendantWeights=np.zeros((2,x.size))
    vref=np.array(vref)
    v_ref.append(np.average(vref[nSteps/2:])*au2wn)
    pop_0.append(np.average(pop[nSteps/2:]))
    v_ref_array_0[n,:]=vref*au2wn
    pop_array_0[n,:]=pop
    PsiHist,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),density=True)
    print n,'from initial propagation of wvf on ground state with pop control v_ref:', v_ref[-1], 'pop:',pop_0[-1]
    #for psi squared
    #first the usual way, adjusting vref during simulation for prop on the ground state
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn0.propagate(x,nDesSteps,nSize=nWalkers)
        v_ref_desc_array_0[0,0,n,:]=np.array(vrd)*au2wn
        pop_desc_array_0[0,0,n,:]=popd
        desc_weight_array_0[i]=descendantWeights

    AveDescendantWeightsGroundPop.append(np.average(desc_weight_array_0,axis=0))
    plt.scatter(x,AveDescendantWeightsGroundPop[n])
    plt.show()
    Psi2Hist_0pop_0pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),
                                              density=True,weights=AveDescendantWeightsGroundPop[n])
                                          
    #second the unusual way, setting vref to be the 'proper' value for prop on the ground state
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn0.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref_0)
        v_ref_desc_array_0[0,1,n,:]=np.array(vrd)*au2wn
        pop_desc_array_0[0,1,n,:]=popd
        desc_weight_array_0[i]=descendantWeights

                                          
    AveDescendantWeightsGroundVRef.append(np.average(desc_weight_array_0,axis=0))
    Psi2Hist_0pop_0vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),
                                               density=True,weights=AveDescendantWeightsGroundVRef[n])

    xhists_0[0,n]=PsiHist
    xhists_0[1,n]=Psi2Hist_0pop_0pop
    xhists_0[2,n]=Psi2Hist_0pop_0vref

    #first the usual way, adjusting vref during simulation for prop on the excited state
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Right.propagate(x,nDesSteps,nSize=nWalkers/2)
        v_ref_desc_array_0[1,0,n,:]=np.array(vrd)*au2wn
        pop_desc_array_0[1,0,n,:]=popd
        desc_weight_array_0[i]=descendantWeights
    AveDescendantWeightsExcRightPop.append(np.average(desc_weight_array_0,axis=0))
    Psi2Hist_0pop_1pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),
                                              density=True,weights=AveDescendantWeightsExcRightPop[n])
    #second the unusual way, setting vref to be the 'proper' value for prop on the excited state
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Right.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref_1)
        v_ref_desc_array_0[1,1,n,:]=np.array(vrd)*au2wn
        pop_desc_array_0[1,1,n,:]=popd
        desc_weight_array_0[i]=descendantWeights
    AveDescendantWeightsExcRightVRef.append(np.average(desc_weight_array_0,axis=0))
    Psi2Hist_0pop_1vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),
                                               density=True,weights=AveDescendantWeightsExcRightVRef[n])
    
    #first the usual way, adjusting vref during simulation for prop on the excited state
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Left.propagate(x,nDesSteps,nSize=nWalkers/2)
        v_ref_desc_array_0[2,0,n,:]=np.array(vrd)*au2wn
        pop_desc_array_0[2,0,n,:]=popd
        desc_weight_array_0[i]=descendantWeights
    AveDescendantWeightsExcLeftPop.append(np.average(desc_weight_array_0,axis=0))
    Psi2Hist_0pop_2pop,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),
                                               density=True,weights=AveDescendantWeightsExcLeftPop[n])   

    #second the unusual way, setting vref to be the 'proper' value for prop on the excited state
    for i in range(nRepsDesc):
        vrd,popd,xd,descendantWeights=Wfn1Left.propagate(x,nDesSteps,setV_ref=True,ConstantV_ref=idealVref_1)
        v_ref_desc_array_0[2,1,n,:]=np.array(vrd)*au2wn
        pop_desc_array_0[2,1,n,:]=popd
        desc_weight_array_0[i]=descendantWeights
    AveDescendantWeightsExcLeftVRef.append(np.average(desc_weight_array_0,axis=0))
    Psi2Hist_0pop_2vref,bin_edges=np.histogram(x, bins=nBins, range=(-1.0,1.0),
                                               density=True,weights=AveDescendantWeightsExcLeftVRef[n])   


    xhists_0[3,n]=Psi2Hist_0pop_1pop-Psi2Hist_0pop_2pop
    xhists_0[4,n]=Psi2Hist_0pop_1vref-Psi2Hist_0pop_2vref
    

    AveDescendantWeightsExcPop.append(AveDescendantWeightsExcRightPop[n]-AveDescendantWeightsExcLeftPop[n])
    AveDescendantWeightsExcVRef.append(AveDescendantWeightsExcRightVRef[n]-AveDescendantWeightsExcLeftVRef[n])
    nonZeros=np.logical_not((AveDescendantWeightsGroundPop[n]==0.0))

    #<0|R2|0>/<0|0> for pop
    numerator=np.sum(x[nonZeros]*x[nonZeros]*
		     ((AveDescendantWeightsGroundPop[n][nonZeros])**2)/
		     AveDescendantWeightsGroundPop[n][nonZeros])
    denominator=np.sum((AveDescendantWeightsGroundPop[n][nonZeros]**2)/
		       AveDescendantWeightsGroundPop[n][nonZeros])
    R2_0Pop[n]=numerator/denominator
    #<1|R2|1>/<1|1> for pop
    numerator=np.sum(x[nonZeros]*x[nonZeros]*
                     ((AveDescendantWeightsExcPop[n][nonZeros])**2)/
		     AveDescendantWeightsGroundPop[n][nonZeros])
    denominator=np.sum((AveDescendantWeightsExcPop[n][nonZeros]**2)/
		       AveDescendantWeightsGroundPop[n][nonZeros])
    R2_1Pop[n]=numerator/denominator

    #<0|R2|0>/<0|0> for vref
    nonZeros=np.logical_not((AveDescendantWeightsGroundVRef[n]==0.0))
    numerator=np.sum(x[nonZeros]*x[nonZeros]*
                     ((AveDescendantWeightsGroundVRef[n][nonZeros])**2)/
                     AveDescendantWeightsGroundVRef[n][nonZeros])
    denominator=np.sum((AveDescendantWeightsGroundVRef[n][nonZeros]**2)/
                       AveDescendantWeightsGroundVRef[n][nonZeros])
    R2_0VRef[n]=numerator/denominator
    #<1|R2|1>/<1|1> for vref
    numerator=np.sum(x[nonZeros]*x[nonZeros]*
                     ((AveDescendantWeightsExcVRef[n][nonZeros])**2)/
                     AveDescendantWeightsGroundVRef[n][nonZeros])
    denominator=np.sum((AveDescendantWeightsExcVRef[n][nonZeros]**2)/
                       AveDescendantWeightsGroundVRef[n][nonZeros])
    R2_1VRef[n]=numerator/denominator


    np.savetxt(Destination+'AverageWeights-n-'+str(n)+'.data',zip(
            AveDescendantWeightsGroundPop[n], 
            AveDescendantWeightsExcRightPop[n],AveDescendantWeightsExcLeftPop[n],AveDescendantWeightsExcPop[n],
            AveDescendantWeightsExcRightVRef[n],AveDescendantWeightsExcLeftVRef[n],AveDescendantWeightsExcVRef[n]))
    
                                           
print 'average energy',np.average(v_ref)
print 'standard deviation', np.std(v_ref)
print 'uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print 'average population', np.average(pop_0)

print 'Analytical<0|R^2|0> is', 0.5/Wfn0.getAlpha()
#print 'Numerical', R2_LeftPop,R2_LeftVRef, R2_RightPop,R2_RightVRef                                                              
print 'averaged', np.average(R2_0Pop), np.average(R2_0VRef)
print 'standard deviation', np.std(R2_0Pop),np.std(R2_0VRef)
unc=[(np.max(R2_0Pop)-np.min(R2_0Pop))/(2.0*np.sqrt(nReps)),(np.max(R2_0VRef)-np.min(R2_0VRef))/(2.0*np.sqrt(nReps))]
print 'uncertainity', unc
fileR2Data=open('R2Data-v0tov1.data','a')
fileR2Data.write(str(0)+'   '+str(nWalkers)+'   '+str(nSteps)+'   '+str(nReps)+'   '+str(nDesSteps)+'   '+str(nRepsDesc)+'      ')
fileR2Data.write(str(0.5*1.0/Wfn0.getAlpha())+'   '+str(np.average(R2_0Pop))+'   '+str(np.average(R2_0VRef))+'   '+str(np.std(R2_0Pop))+'   '+str(np.std(R2_0VRef))+'   '+str(unc[0])+'   '+str(unc[1])+'\n')



print 'Analytical<1|R^2|1> is', 1.5/Wfn0.getAlpha()

print 'averaged', np.average(R2_1Pop), np.average(R2_1VRef)
print 'standard deviation', np.std(R2_1Pop),np.std(R2_1VRef)
unc=[(np.max(R2_1Pop)-np.min(R2_1Pop))/(2.0*np.sqrt(nReps)),(np.max(R2_1VRef)-np.min(R2_1VRef))/(2.0*np.sqrt(nReps))]
print 'uncertainity', unc
fileR2Data=open('R2Data-v0tov1.data','a')
fileR2Data.write(str(1)+'   '+str(nWalkers)+'   '+str(nSteps)+'   '+str(nReps)+'   '+str(nDesSteps)+'   '+str(nRepsDesc)+'      ')
fileR2Data.write(str(1.5*1.0/Wfn0.getAlpha())+'   '+str(np.average(R2_1Pop))+'   '+str(np.average(R2_1VRef))+'   '+str(np.std(R2_1Pop))+'   '+str(np.std(R2_1VRef))+'   '+str(unc[0])+'   '+str(unc[1])+'\n')


fileR2Data.close()


np.savetxt(Destination+'vrefForPOP.data', v_ref_array_0.T)
np.savetxt(Destination+'popForPOP.data', pop_array_0.T)
np.savetxt(Destination+'vrefDescFor-POP0-POP0.data',v_ref_desc_array_0[0,0].T)
np.savetxt(Destination+'vrefDescFor-POP0-VREF0.data',v_ref_desc_array_0[0,1].T)
np.savetxt(Destination+'popDescFor-POP0-VREF0.data',pop_desc_array_0[0,1].T)
np.savetxt(Destination+'popDescFor-POP0-POP0.data',pop_desc_array_0[0,0].T)
np.savetxt(Destination+'vrefDescFor-POP0-POP1.data',v_ref_desc_array_0[1,0].T)
np.savetxt(Destination+'vrefDescFor-POP0-VREF1.data',v_ref_desc_array_0[1,1].T)
np.savetxt(Destination+'popDescFor-POP0-VREF1.data',pop_desc_array_0[1,1].T)
np.savetxt(Destination+'popDescFor-POP0-POP1.data',pop_desc_array_0[1,0].T)
np.savetxt(Destination+'vrefDescFor-POP0-POP2.data',v_ref_desc_array_0[2,0].T)
np.savetxt(Destination+'vrefDescFor-POP0-VREF2.data',v_ref_desc_array_0[2,1].T)
np.savetxt(Destination+'popDescFor-POP0-VREF2.data',pop_desc_array_0[2,1].T)
np.savetxt(Destination+'popDescFor-POP0-POP2.data',pop_desc_array_0[2,0].T)

print 'ideal is',idealVref_1*au2wn
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

colors=['r','b','k','g']
#plotting first v_Ref

for n in range(nReps):
    plt.scatter(xCoord0[n],AveDescendantWeightsGroundPop[n],c=red)
    plt.scatter(xCoord0[n],AveDescendantWeightsGroundVRef[n],c=blue)
plt.ylabel('Average Descendant Weight d_0')
plt.xlabel('Coordinate x (bohr)')
plt.savefig(Destination+'UWplotGroundState.png')
plt.clf()

for n in range(nReps):
    plt.scatter(xCoord0[n],AveDescendantWeightsExcPop[n],c=red)
    plt.scatter(xCoord0[n],AveDescendantWeightsExcVRef[n],c=blue)
plt.ylabel('Average Descendant Weight d_1')
plt.xlabel('Coordinate x (bohr)')
plt.savefig(Destination+'UWplotExcitedState.png')
plt.clf()

for n in range(nReps):
    plt.plot(np.arange(nSteps),v_ref_array_0[n],c=red)
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[0,0,n],c=orange)
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[0,1,n],c=darkOrange)
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[1,0,n],c=purple)
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[1,1,n],c=darkPurple)
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[2,0,n],c=pink)
    plt.plot(np.arange(nDesSteps)+nSteps,v_ref_desc_array_0[2,1,n],c=darkPink)
plt.ylabel('V_Ref (1/cm)')
plt.ylim(0,6000)
plt.xlabel('Simulation Time Step (5 a.u.)')
plt.savefig(Destination+'vRef.png')
plt.clf()    
for n  in range(nReps):
    plt.plot(np.arange(nSteps),pop_array_0[n],c=red)
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[0,0,n],c=orange)
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[0,1,n],c=darkOrange)
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[1,0,n],c=purple)
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[1,1,n],c=darkPurple)
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[2,0,n],c=pink)
    plt.plot(np.arange(nDesSteps)+nSteps,pop_desc_array_0[2,1,n],c=darkPink)
plt.ylabel('Population of Walkers')
plt.xlabel('Simulation Time Step ('+str(Wfn0.get_dtau())+' a.u.)')
plt.savefig(Destination+'population.png')
plt.clf()

bin_center=(bin_edges[:-1]+bin_edges[1:])/2.0
for n in range(nReps):
    plt.plot(bin_center,xhists_0[0,n],c=red,linewidth=.5)
plt.plot(bin_center,np.average(xhists_0[0],axis=0),c=darkRed,linewidth=2.0)
plt.ylabel('Psi(x)')
plt.xlabel('x (bohr)')
plt.savefig(Destination+"Psis.png")
plt.clf()

for n in range(nReps):
    plt.plot(bin_center,xhists_0[1,n],c=red,linewidth=.1)
    plt.plot(bin_center,xhists_0[2,n],c=darkRed,linewidth=.1)
    plt.plot(bin_center,xhists_0[3,n],c=purple,linewidth=.1)
    plt.plot(bin_center,xhists_0[4,n],c=darkPurple,linewidth=.1)
plt.plot(bin_center,np.average(xhists_0[1],axis=0),c=red,linewidth=2)
plt.plot(bin_center,np.average(xhists_0[2],axis=0),c=darkRed,linewidth=2)
plt.plot(bin_center,np.average(xhists_0[3],axis=0),c=purple,linewidth=2)
plt.plot(bin_center,np.average(xhists_0[4],axis=0),c=darkPurple,linewidth=2)
plt.ylabel('Psi_n(x)*Psi_m(x)')
plt.xlabel('x (bohr)')
plt.savefig(Destination+'Psi-star-psi.png')
plt.clf()

np.savetxt(Destination+'PsiHist-pop.data',xhists_0[0].T)
np.savetxt(Destination+'PsiHist-AVERAGE-pop.data',zip(bin_center,np.average(xhists_0[0],axis=0)))
np.savetxt(Destination+'Psi2Hist-pop-pop.data',xhists_0[1].T)
np.savetxt(Destination+'Psi2Hist-AVERAGE-pop-pop.data',zip(bin_center,np.average(xhists_0[1],axis=0)))
np.savetxt(Destination+'Psi2Hist-pop-vref.data',xhists_0[2].T)
np.savetxt(Destination+'Psi2Hist-AVERAGE-pop-vref.data',zip(bin_center,np.average(xhists_0[2],axis=0)))

print 'done!'

