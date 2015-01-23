#include <iomanip>
#include <iostream>
#include <string>
#include "BranchChanges.h"
#include "CharacterChanges.h"
#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "Settings.h"
#include "Tree.h"


Mcmc::Mcmc(Settings* sp, Model* mp, MbRandom* rp) {

    modelPtr        = mp;
    ranPtr          = rp;
    settingsPtr     = sp;
    
    chainLength     = settingsPtr->getChainLength();
    burnIn          = settingsPtr->getBurnIn();
    printFrequency  = settingsPtr->getPrintFrequency();
    sampleFrequency = settingsPtr->getSampleFrequency();

	std::string parmFileName = settingsPtr->getOutputFileName() + ".parm";
	std::string treeFileName = settingsPtr->getOutputFileName() + ".tree";
	parmOut.open( parmFileName.c_str(), std::ios::out );
	treeOut.open( treeFileName.c_str(), std::ios::out );
    
    runChain();
    
    parmOut.close();
    treeOut.close();
}

Mcmc::~Mcmc(void) {

}

void Mcmc::runChain(void) {

    int curSpace = 0;
    int newSpace = 1;
    double curLnL = modelPtr->lnLikelihood(curSpace);
    for (int n=1; n<=chainLength; n++)
        {
        // propose a new value for the Markov chain
        double lnProposalRatio = modelPtr->update(newSpace);
        
        // calculate the acceptance probability
        double newLnL = modelPtr->lnLikelihood(newSpace);
        double lnLikelihoodRatio = newLnL - curLnL;
        double lnPriorRatio = modelPtr->lnPriorRatio(newSpace, curSpace);
        double lnR = lnLikelihoodRatio + lnPriorRatio + lnProposalRatio;

        if (n % printFrequency == 0)
            {
            std::cout << std::fixed << std::setprecision(8);
            std::cout << n << " -- " << curLnL << " -> " << newLnL;
            }
            
        // accept or reject
        if ( log(ranPtr->uniformRv()) < lnR )
            {
            // accept
            modelPtr->accept();
            curLnL = newLnL;
            if (n % printFrequency == 0)
                std::cout << " (Accept) " << modelPtr->getPi(0)[0] << " " << modelPtr->averageBranchLength(0) << std::endl;
            }
        else
            {
            // reject
            modelPtr->reject();
            if (n % printFrequency == 0)
                std::cout << " (Reject) " << modelPtr->getPi(0)[0] << " " << modelPtr->averageBranchLength(0) << std::endl;
            }
            
        // sample chain
        if ( n >= burnIn && n % sampleFrequency == 0 )
            modelPtr->mapCharacters(0);
        if ( n % sampleFrequency == 0 || n == 1 || n == chainLength )
            saveStates(n, curLnL);
        }
}

void Mcmc::saveStates(int n, double lnL) {

    // save to the parameter file
    if (n == 1)
        parmOut << "Cycle" << '\t' << "lnL" << '\t' << "pi[0]" << '\t' << "pi[1]" << '\t' << "Sum(V)" << '\n';
    double* pi = modelPtr->getPi(0);
    parmOut << n << '\t' << lnL << '\t' << pi[0] << '\t' << pi[1] << '\t' << modelPtr->averageBranchLength(0) << '\n';

    // save to tree file
    std::string treeStr = modelPtr->getTree(0)->getTreeString();
    treeOut << treeStr << '\n';
    
    // save to character mapping file (this file is over-written each time)
    if ( n % (sampleFrequency*10) == 0 || n == chainLength )
        {
        std::ofstream charOut;
        std::string charFileName = settingsPtr->getOutputFileName() + ".char";
        charOut.open( charFileName.c_str(), std::ios::out );
        CharacterChanges* c0 = modelPtr->getCharacterChange(0);
        std::string header = c0->getChangeProbsHeader();
        charOut << "Character" << '\t' << header << '\n';
        for (int i=0; i<modelPtr->getNumGenes(); i++)
            {
            CharacterChanges* ci = modelPtr->getCharacterChange(i);
            std::string dataStr = ci->getChangeProbs();
            charOut << i+1 << '\t' << dataStr << '\n';
            }
        charOut.close();
        }

    // also write files with the information for each node
    if ( n % (sampleFrequency*10) == 0 || n == chainLength )
        {
        for (int n=0; n<modelPtr->getTree(0)->getNumNodes(); n++)
            {
            BranchChanges* bc = modelPtr->getBranchChangesForNode(n);
            
            std::string nodeName = bc->getName();
            char cStr[100];
            if (nodeName == "")
                sprintf(cStr, "node_%d", n);
            else
                sprintf(cStr, "node_%s", nodeName.c_str());
            std::string branchFileName = settingsPtr->getOutputFileName() + "." + cStr;

            std::ofstream branchOut;
            branchOut.open( branchFileName.c_str(), std::ios::out );
            for (int i=0; i<modelPtr->getNumGenes(); i++)
                {
                branchOut << i+1 << '\t';
                branchOut << bc->getAverageNumLosses(i) << '\t';
                branchOut << bc->getAverageNumGains(i) << '\t';
                branchOut << bc->getProbabilityOne(i) << '\n';
                }
            branchOut.close();
            }
        }
}


