'''
Randomly sample from samples to identifiy how many networks are needed
to build robust communities
'''


import panCancerTest as pct

#looad cancer data
patDiffs = pct.loadCancerData(0.01)


def buildHyphaFromSampledData(allDat, numSamps):
    #this builds a hypha from a sampling of data
    this_hyp = None
    return this_hyp


#sample patient data from 5 to all for 100 times each
#compute NMI between hyphae and original graph
#compute NMI between cancer types

#compute mean NMI
