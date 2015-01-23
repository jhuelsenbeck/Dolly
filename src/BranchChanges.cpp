#include "BranchChanges.h"


BranchChanges::BranchChanges(int nc, int idx, std::string nme) {

    numCharacters = nc;
    myNodeIndex   = idx;
    name          = nme;
    
    numSamples   = new int[numCharacters];
    numGains     = new int[numCharacters];
    numLosses    = new int[numCharacters];
    numOneAtNode = new int[numCharacters];
    
    for (int i=0; i<numCharacters; i++)
        {
        numSamples[i]   = 0;
        numGains[i]     = 0;
        numLosses[i]    = 0;
        numOneAtNode[i] = 0;
        }
}

BranchChanges::~BranchChanges(void) {

    delete [] numSamples;
    delete [] numGains;
    delete [] numLosses;
    delete [] numOneAtNode;
}

void BranchChanges::updateChangesForSite(int site, int n01, int n10, int endState) {

    numSamples[site]   += 1;
    numGains[site]     += n01;
    numLosses[site]    += n10;
    numOneAtNode[site] += endState;
}

double BranchChanges::getAverageNumGains(int idx) {

    return (double)numGains[idx] / numSamples[idx];
}

double BranchChanges::getAverageNumLosses(int idx) {

    return (double)numLosses[idx] / numSamples[idx];
}

double BranchChanges::getProbabilityOne(int idx) {

    return (double)numOneAtNode[idx] / numSamples[idx];
}


