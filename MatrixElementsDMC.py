import  dmc1D as dmc
import numpy as np
nWalkers=2000
Wfn0=dmc.wavefunction(nWalkers,'harmonic',plotting=False)
v_ref=[]
nReps=10
for n in range(nReps):
    v_ref.append(Wfn0.propagate(Wfn0.xcoords,100))
print 'average energy',np.average(v_ref)
print 'standard deviation', np.std(v_ref)
print 'uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print Wfn0.omega/2.0

Wfn1=dmc.wavefunction(nWalkers,'half harmonic',plotting=False)
v_ref_1=[]
for n in range(nReps):
    v_ref.append(Wfn1.propagate(Wfn1.xcoords,100))
print 'v1 average energy',np.average(v_ref)
print 'v1 standard deviation', np.std(v_ref)
print 'v1 uncertainity', (np.max(v_ref)-np.min(v_ref))/(2.0*np.sqrt(nReps))
print 3.0*Wfn1.omega/2.0


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



