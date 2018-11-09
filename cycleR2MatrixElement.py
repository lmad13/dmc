import os
import time
nWalkers=2000
nReps=10
nSteps=100
nDescSteps=25
nRepsDesc=25
count=0
starttime=time.time()
timeRate=4.8E-6  #my rough estimate
for nDescSteps in [10,100,200,500]:
    for nRepsDesc in [1,2]:
        for nReps in [1,2]:
            for nWalkers in [2000,2000]:
                count=count+1
                print '\n\n################    WORKING ON #',count,
                print str(len(range(50,125,25))*len(range(25,125,25))*len((5,15,5))*len(range(1000,1500,500))),' #########'
                print 'nWalkers:   ',nWalkers
                print 'nReps:      ',nReps
                print 'nRepsDesc:  ',nRepsDesc
                print 'nDescSteps:  ',nDescSteps
                print 'nSteps:     ',nSteps
                print '#####   Estimated Time:',str(timeRate*nWalkers*nReps*nDescSteps*nRepsDesc*nSteps),' s #######\n'
                executeStartTime=time.time()
                os.system('./R2MatrixElementv0tov1.py '+str(nWalkers)+' '+str(nReps)+
                          ' '+str(nSteps)+' '+str(nDescSteps)+' '+str(nRepsDesc))
                executionTime=time.time()-executeStartTime
                timeRate=executionTime/(nWalkers*nReps*nDescSteps*nRepsDesc*nSteps)
                print 'that took', executionTime, 'with a timerate of ', timeRate

                
                
print 'done cycling!', time.time()-starttime
