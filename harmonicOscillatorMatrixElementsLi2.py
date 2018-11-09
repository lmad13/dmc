import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys
au2wn=219474.63
nWalkers=20000

Wfn0=dmc.wavefunction(nWalkers,'harmonic',plotting=False,omegaInput=913.65,molecule='Li2')
Wfn1=dmc.wavefunction(nWalkers,'half harmonic',plotting=False,omegaInput=913.65,molecule='Li2')

Wfn0.getIdealVref()
Wfn1.getIdealVref()

print 'the reduced mass is ', Wfn0.mass/1822.88839, 'amu?'

vref_list_0=[]
vref_list_1=[]
pop_0=[]
pop_1=[]
nReps=10
nSteps=2000
nEquilibrationSteps=12000
nDesSteps=75
nRepsDW=2000
print 'important parameters:'
print 'nEquilibrationSteps',nEquilibrationSteps,'\n nSteps',nSteps, 'dtau', Wfn0.dtau
print 'nRepsDW',nRepsDW, '\n nDesSteps',nDesSteps
vref_array_0=np.zeros((nReps,nSteps))
vref_desc_array_0=np.zeros((2,nReps,nDesSteps))
pop_array_0=np.zeros((nReps,nSteps))
pop_desc_array_0=np.zeros((2,nReps,nDesSteps))

vref_array_1=np.zeros((nReps,nSteps))
vref_desc_array_1=np.zeros((2,nReps,nDesSteps))
pop_array_1=np.zeros((nReps,nSteps))
pop_desc_array_1=np.zeros((2,nReps,nDesSteps))

nBins=101

xhists_0=np.zeros((3,nReps,nBins)) #Psi, Psi_popPsi_pop, Psi_popPsi_vref

Wfn0.setX(Wfn0.xcoords+.2)
Wfn1.xcoords[:nWalkers/2]=.2
Wfn1.xcoords[nWalkers/2:]=-.2
if 'new' in sys.argv:
    v,pop,x0Equilibration,d=Wfn0.propagate(Wfn0.xcoords,nEquilibrationSteps)
    v,pop,x1Equilibration,d=Wfn1.propagate(Wfn1.xcoords,nEquilibrationSteps)
    np.savetxt('GroundState20000Equilibration.coord',x0Equilibration)
    np.savetxt('ExcState20000Equilibration.coord',x1Equilibration)
else:
    x0Equilibration=np.loadtxt('GroundState20000Equilibration.coord')
    x1Equilibration=np.loadtxt('ExcState20000Equilibration.coord')

print 'Analytical <0|R^2|0> is', 0.5/Wfn0.getAlpha()
print 'Analytical <1|R^2|1> is ',1.5*1.0/Wfn0.getAlpha()
print 'Analytical <1|R|1> is   ',np.sqrt(1.0/(2.0*Wfn0.getAlpha()))
idealVref0=Wfn0.getIdealVref()
idealVref1=Wfn1.getIdealVref()

print 'The ideal energies', idealVref0*au2wn,idealVref1*au2wn

print '\n\n'
print 'n     E_0         <0|x^2|0>        E_1      <1|x^2|1>       <<0>|x^2|0>       <0|x|1>'
for n in range(nReps):
    #Equilibrate on Ground
    vref_0,pop,x0,d=Wfn0.propagate(x0Equilibration,nSteps)
    vref_0=np.array(vref_0)
    vref_list_0.append(np.average(vref_0[nSteps/2:])*au2wn)
    pop_0.append(np.average(pop[nSteps/2:]))
    vref_array_0[n,:]=vref_0*au2wn
    pop_array_0[n,:]=pop
    
    vrd,popd,xd,descendantWeights0=Wfn0.propagate(x0,nDesSteps,nSize=nWalkers)

    vref_desc_array_0[0,n,:]=np.array(vrd)*au2wn


    #Equilibrate on the Excited State
    vref_1,pop,x1,d=Wfn1.propagate(x1Equilibration,nSteps)
    vref_1=np.array(vref_1)
    vref_list_1.append(np.average(vref_1[nSteps/2:])*au2wn)
    pop_1.append(np.average(pop[nSteps/2:]))
    vref_array_1[n,:]=vref_1*au2wn
    
    pop_array_1[n,:]=pop
    vrd,popd,xd,descendantWeights1=Wfn1.propagate(x1,nDesSteps,nSize=nWalkers)
    print n ,np.average(vref_list_0),
    print  np.sum(x0*x0*descendantWeights0)/np.sum(descendantWeights0),
    print np.average(vref_list_1),
    print  np.sum(x1*x1*descendantWeights1)/np.sum(descendantWeights1),


    #Determine descendant weights on excited and ground states from ground state population via the averaging 
    descendantWeights00=np.zeros((x0.shape[0]))
    descendantWeights01=np.zeros((x0.shape[0]))
    for irep in range(nRepsDW):
         v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn0.propagate(x0,nDesSteps,nSize=nWalkers)
         descendantWeights00=descendantWeights00+descendantsTemp
         v_ref_DW_list,pop_DW_list,DWFinalCoords,descendantsTemp=Wfn1.propagate(x0,nDesSteps,nSize=nWalkers)
         descendantWeights01=descendantWeights01+descendantsTemp


    descendantWeights00=descendantWeights00/nRepsDW
    descendantWeights01=descendantWeights01/nRepsDW

    nonZeros_0=np.logical_not((descendantWeights00==0.0))
    print np.sum(x0[nonZeros_0]*x0[nonZeros_0]*descendantWeights00[nonZeros_0])/np.sum(descendantWeights00[nonZeros_0]),
    print  np.sum(x0[nonZeros_0]*x0[nonZeros_0]*(descendantWeights01[nonZeros_0]*descendantWeights01[nonZeros_0])/descendantWeights00[nonZeros_0])/np.sum((descendantWeights01[nonZeros_0]*descendantWeights01[nonZeros_0])/descendantWeights00[nonZeros_0]),
    
    print  np.sum(np.abs(x0[nonZeros_0])*descendantWeights01[nonZeros_0])/np.sqrt(np.sum(descendantWeights00[nonZeros_0])*np.sum(descendantWeights01[nonZeros_0]**2/descendantWeights00[nonZeros_0]))
    #print np.sum(descendantWeights00[nonZeros_0]), np.sum(descendantWeights01[nonZeros_0]), np.sum(descendantWeights01[nonZeros_0]*descendantWeights01[nonZeros_0]/descendantWeights00[nonZeros_0])
    sys.stdout.flush()
    
    
print 'shape', vref_array_0.shape
print '.....average vref_0', np.average(np.average(vref_array_0,axis=1)),'+/-',np.std(np.average(vref_array_0,axis=1))
print '.....average pop',    np.average(np.average(pop_array_0,axis=1)),'+/-',np.std(np.average(pop_array_0,axis=1))
print '.....average vref_1', np.average(np.average(vref_array_1,axis=1)),'+/-',np.std(np.average(vref_array_1,axis=1))
print '.....average pop',    np.average(np.average(pop_array_1,axis=1)),'+/-',np.std(np.average(pop_array_1,axis=1))

print 'Analytical <0|R^2|0> is', 0.5/Wfn0.getAlpha()
print 'Analytical <1|R^2|1> is ',1.5*1.0/Wfn0.getAlpha()
print 'Analytical <1|R|1> is   ',np.sqrt(1.0/(2.0*Wfn0.getAlpha()))
idealVref0=Wfn0.getIdealVref()
idealVref1=Wfn1.getIdealVref()

print 'The ideal energies', idealVref0*au2wn,idealVref1*au2wn

print 'done!'
