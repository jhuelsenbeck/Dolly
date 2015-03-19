#include <istream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "AncestralProbabilities.h"
#include "BranchChanges.h"
#include "CharacterChanges.h"
#include "Clades.h"
#include "Data.h"
#include "MbRandom.h"
#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "Settings.h"
#include "Tree.h"

#define BAD_LIKE -10e50



Model::Model(Settings* sp, Data* dp, MbRandom* rp) {
    
    std::cout << "   * Setting up the phylogenetic model" << std::endl;

    // remember important objects
    dataPtr     = dp;
    settingsPtr = sp;
    ranPtr      = rp;
    
    // instantiate model parameters
    branchLengthPriorParm = settingsPtr->getBranchLengthPrior();
    alpha.push_back(1.0);
    alpha.push_back(1.0);
    for (int i=0; i<2; i++)
        pi[0][i] = pi[1][i] = 0.5;
    
    theTree[0]  = NULL;
    theTree[1]  = NULL;
    readTree();
    theTree[1] = new Tree(*theTree[0], 1);
    theTree[0]->setMySpace(0);
    theTree[1]->setMySpace(1);
    //theTree[0]->print();
    //theTree[1]->print();
    
    // initialize conditional likelihoods
    initializeConditionalLikelihood();
    
    // set the proposal probabilities
    proposalRates.push_back(20);    // change tree
    proposalRates.push_back(1);     // change pi
    double sum = 0.0;
    for (int i=0; i<proposalRates.size(); i++)
        sum += proposalRates[i];
    for (int i=0; i<proposalRates.size(); i++)
        proposalProbs.push_back(proposalRates[i]/sum);
    
    // allocate information to remember changes
    for (int i=0; i<dataPtr->getNumChar(); i++)
        characterChanges.push_back(new CharacterChanges(i));
    for (int i=0; i<theTree[0]->getNumNodes(); i++)
        {
        Node* p = theTree[0]->getNodeWithIndex(i);
        branchChanges.push_back(new BranchChanges(dataPtr->getNumChar(), i, p->getName()) );
        }
    
    std::cout << "lnL(0) = " << lnLikelihood(0) << std::endl;
    std::cout << "lnL(1) = " << lnLikelihood(1) << std::endl;
}

Model::~Model(void) {

    if (theTree[0] != NULL)
        delete theTree[0];
    if (theTree[1] != NULL)
        delete theTree[1];
    for (int i=0; i<characterChanges.size(); i++)
        delete characterChanges[i];
    for (int i=0; i<branchChanges.size(); i++)
        delete branchChanges[i];
}

void Model::accept(void) {

    *theTree[0] = *theTree[1];
    pi[0][0] = pi[1][0];
    pi[0][1] = pi[1][1];
}

double Model::averageBranchLength(int space) {

    Tree* t = theTree[space];
    int count = 0;
    double sum = 0.0;
    for (int n=0; n<t->getNumberOfDownPassNodes(); n++)
        {
        Node* p = t->getDownPassNode(n);
        if (p != t->getRoot())
            {
            sum += p->getBranchLength();
            count++;
            }
        }
    return sum / count;
}

void Model::calcConditionalLikelihoodsDown(int space) {

    /*
    Tree* t = theTree[space];
    
    for (int n=0; n<t->getNumberOfDownPassNodes(); n++)
        {
        Node* p = t->getDownPassNode(n);
        if (p->isLeaf() == false)
            {
            std::set<Node*> descendants = p->getDescendants();
            int numUnflaggedDescendants = p->numberOfDescendantsWithFlag(false);
            double* clP = &cls[space][p->getIndex()][0];
            std::vector<double*> clQ(numUnflaggedDescendants);
            std::vector<double**> ti(numUnflaggedDescendants);
            int k = 0;
            for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
                {
                if ( (*it)->getFlag() == false )
                    {
                    clQ[k] = &cls[space][(*it)->getIndex()][0];
                    ti[k] = (*it)->exposeTransitionProbability();
                    k++;
                    }
                }
                
            for (int c=0; c<numGenes; c++)
                {
                for (int i=0; i<2; i++)
                    {
                    double product = 1.0;
                    for (int d=0; d<clQ.size(); d++)
                        {
                        double sum = 0.0;
                        for (int j=0; j<2; j++)
                            sum += ti[d][i][j] * clQ[d][j];
                        product *= sum;
                        }
                    clP[i] = product;
                    }
                    
#               if 0
                if (c == 0 || c == 100)
                    std::cout << "D " << c << " " << p->getIndex() << " " << clP[0] << " " << clP[1] << std::endl;
#               endif
                    
                clP += 2;
                for (int d=0; d<clQ.size(); d++)
                    clQ[d] += 2;
                }
            }
        }*/
}

void Model::calcConditionalLikelihoodsUpToNode(int space, Node* theNode) {
    
    /*
    Tree* t = theTree[space];
    
    for (int n=t->getNumberOfDownPassNodes()-1; n>=0; n--)
        {
        Node* p = t->getDownPassNode(n);
        if (p->isLeaf() == false && p->getAncestor() != NULL)
            {
            std::set<Node*> descendants = p->getDescendants();
            int numUnflaggedDescendants = p->numberOfDescendantsWithFlag(false);
            numUnflaggedDescendants++;
            double* clP = &cls[space][p->getIndex()][0];
            std::vector<double*> clQ(numUnflaggedDescendants);
            std::vector<double**> ti(numUnflaggedDescendants);
            int k = 0;
            for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
                {
                if ( (*it)->getFlag() == false )
                    {
                    clQ[k] = &cls[space][(*it)->getIndex()][0];
                    ti[k] = (*it)->exposeTransitionProbability();
                    k++;
                    }
                }
            clQ[k] = &cls[space][p->getAncestor()->getIndex()][0];
            ti[k] = p->exposeTransitionProbability();
                
            for (int c=0; c<numGenes; c++)
                {
                for (int i=0; i<2; i++)
                    {
                    double product = 1.0;
                    for (int d=0; d<clQ.size(); d++)
                        {
                        double sum = 0.0;
                        for (int j=0; j<2; j++)
                            sum += ti[d][i][j] * clQ[d][j];
                        product *= sum;
                        }
                    clP[i] = product;
                    }

#               if 0
                if (c == 0 || c == 100)
                    {
                    std::cout << "U " << std::setw(5) << c << " " << p->getIndex() << " " << clP[0] << " " << clP[1] << std::endl;
                    for (int d=0; d<clQ.size(); d++)
                        std::cout << "       " << clQ[d][0] << " " << clQ[d][1] << std::endl;
                    }
#               endif

                clP += 2;
                for (int d=0; d<clQ.size(); d++)
                    clQ[d] += 2;
                }
            }
        }*/
}

void Model::calculateStateProbabilitiesAtNode(int space, Node* theNode) {

#   if 0
    double* clP = &cls[space][theNode->getIndex()][0];

    for (int c=0; c<numGenes; c++)
        {
        double sum = 0.0;
        for (int i=0; i<2; i++)
            sum += clP[i] * pi[i];
            
        for (int i=0; i<2; i++)
            clP[i] *= (pi[i] / sum);

#       if 1
        if (1 /*c == 0 || c == 100*/)
            {
            std::cout << theNode << " " << c << " -- ";
            for (int i=0; i<2; i++)
                std::cout << clP[i] << " ";
            std::vector<int> x = dataPtr->getSitePattern(c);
            for (int z=0; z<x.size(); z++)
                std::cout << x[z];
            std::cout << std::endl;
            }
#       endif

        clP += 2;
        }
#   endif
}

void Model::calculateStateProbs(int space) {

    if ( nodesOfInterest.size() == 0 )
        Msg::error("No nodes of interest");
    
    for (std::vector<Node*>::iterator it = nodesOfInterest.begin(); it != nodesOfInterest.end(); it++)
        {
        theTree[space]->markPathToRootFromNode(*it);
        calcConditionalLikelihoodsDown(space);
        calcConditionalLikelihoodsUpToNode(space, *it);
        calculateStateProbabilitiesAtNode(space, *it);
        }

    /*for (int n=0; n<theTree->getNumberOfDownPassNodes(); n++)
        {
        Node* p = theTree->getDownPassNode(n);
        AncestralProbabilities* ap = new AncestralProbabilities(p, numGenes);
        ancestralProbabilities.push_back(ap);
        }*/
}

void Model::initializeConditionalLikelihood(void) {

    numGenes  = dataPtr->getNumChar();
    numNodes  = theTree[0]->getNumNodes();
    numStates = 2;
    //std::cout << "numGenes=" << numGenes << " numNodes=" << numNodes << std::endl;
    
    for (int s=0; s<2; s++)
        {
        std::vector<Node*> treeNodes = theTree[s]->exposeNodes();
        for (int i=0; i<dataPtr->getNumTaxa(); i++)
            {
            Node* nde = treeNodes[i];
            double* p = nde->exposeConditionalLikelihoods();
            for (int j=0; j<numGenes; j++)
                {
                int charCode = dataPtr->getCharacter(i, j);
                int possibleChars[2] = { 0, 0 };
                dataPtr->getPossibleChars(charCode, possibleChars);
                p[0] = (double)possibleChars[0];
                p[1] = (double)possibleChars[1];
                if (settingsPtr->getUseGenePresenceProbs() == true)
                    {
                    if (possibleChars[0] == 0 && possibleChars[1] == 1)
                        {
                        double e = dataPtr->getGenePresenceProbability(i, j);
                        if (e < 0.0 || e > 1.0)
                            Msg::error("Nonsense gene presence probability");
                        p[0] = 1.0 - e;
                        p[1] = e;
                        }
                    }
		//std::cout << "Model::initializeConditionalLikelihood\ttaxon\t" << i << "\tgene\t" << j << "\tp[0]\t" <<
		//  p[0] << "\tp[1]\t" << p[1] << std::endl;
                p += 2;
                }

	    // update transition probabilities if presence probabilities have been incorporated
	    if (settingsPtr->getUseGenePresenceProbs() == true && nde != theTree[s]->getRoot())
	        {
		double x = nde->getBranchProportion();
		nde->setBranchProportion(x);
		}
                
            double* dp = nde->exposeDummyConditionalLikelihoods();
            dp[0] = 1.0;
            dp[1] = 0.0;
            dp[2] = 0.0;
            dp[3] = 1.0;
            }
        }
    
#   if 0
    for (int n=0; n<dataPtr->getNumTaxa(); n++)
        {
        Node* nde0 = theTree[0]->exposeNodes()[n];
        Node* nde1 = theTree[1]->exposeNodes()[n];
        double* p0 = nde0->exposeConditionalLikelihoods();
        double* p1 = nde1->exposeConditionalLikelihoods();

        std::cout << std::setw(5) << n << " -- ";
        for (int i=0; i<numGenes; i++)
            {
            std::cout << std::fixed << std::setprecision(1) << p0[0] << "," << p0[1] << " ";
            p0 += 2;
            }
        std::cout << std::endl << "         ";
        for (int i=0; i<numGenes; i++)
            {
            std::cout << std::fixed << std::setprecision(1) << p1[0] << "," << p1[1] << " ";
            p1 += 2;
            }
        std::cout << std::endl;
        }
#   endif
}

double Model::lnLikelihood(void) {

    return lnLikelihood(0);
}

double Model::lnLikelihood(int space) {

    //return 0.0;
    
    if (pi[space][0] < 0.001 || pi[space][0] > 0.999)
        return BAD_LIKE;
    
    // get a pointer to the tree
    Tree* t = theTree[space];
    
    // pass down the tree initializing the conditional likelihoods
    for (int n=0; n<t->getNumberOfDownPassNodes(); n++)
        {
        Node* p = t->getDownPassNode(n);
        if (p->isLeaf() == false)
            {
            std::set<Node*> descendants = p->getDescendants();
            double* clP = p->exposeConditionalLikelihoods();
            std::vector<double*> clQ(descendants.size());
            std::vector<double**> ti(descendants.size());
            int k = 0;
            for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
                {
                clQ[k] = (*it)->exposeConditionalLikelihoods();
                ti[k] = (*it)->exposeTransitionProbability();
                k++;
                }
                
            for (int c=0; c<numGenes; c++)
                {
                for (int i=0; i<2; i++)
                    {
                    double product = 1.0;
                    for (int d=0; d<clQ.size(); d++)
                        {
                        double sum = 0.0;
                        for (int j=0; j<2; j++)
                            sum += ti[d][i][j] * clQ[d][j];
                        product *= sum;
                        }
                    clP[i] = product;
                    }
                clP += 2;
                for (int d=0; d<clQ.size(); d++)
                    clQ[d] += 2;
                }
                
            }
        }


    // calculate the likelihood
    double lnProbVariable = log(1.0 - probInvariant(space));
    Node* r = t->getRoot();
    double* clP = r->exposeConditionalLikelihoods();
    double lnL = 0.0;
    for (int c=0; c<numGenes; c++)
        {
        double like = pi[space][0] * clP[0] + pi[space][1] * clP[1];
        lnL += log(like) - lnProbVariable;
        clP += 2;
        }
    
    return lnL;
}

double Model::probInvariant(int space) {

    //return 0.0;
#   if defined(CONDITION_NO_CONSTANT)
    int numIllegalPats = 2;
#   elif defined(CONDITION_NO_ZEROS)
    int numIllegalPats = 1;
#   endif

    // get a pointer to the tree
    Tree* t = theTree[space];
    
    // pass down the tree initializing the conditional likelihoods
    for (int n=0; n<t->getNumberOfDownPassNodes(); n++)
        {
        Node* p = t->getDownPassNode(n);
        if (p->isLeaf() == false)
            {
            std::set<Node*> descendants = p->getDescendants();
            double* clP = p->exposeDummyConditionalLikelihoods();
            std::vector<double*> clQ(descendants.size());
            std::vector<double**> ti(descendants.size());
            int k = 0;
            for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
                {
                clQ[k] = (*it)->exposeDummyConditionalLikelihoods();
                ti[k] = (*it)->exposeTransitionProbability();
                k++;
                }
                
            for (int c=0; c<numIllegalPats; c++)
                {
                for (int i=0; i<2; i++)
                    {
                    double product = 1.0;
                    for (int d=0; d<clQ.size(); d++)
                        {
                        double sum = 0.0;
                        for (int j=0; j<2; j++)
                            sum += ti[d][i][j] * clQ[d][j];
                        product *= sum;
                        }
                    clP[i] = product;
                    }
                clP += 2;
                for (int d=0; d<clQ.size(); d++)
                    clQ[d] += 2;
                }
                
            }
        }

    // calculate the likelihood
    Node* r = t->getRoot();
    double* clP = r->exposeDummyConditionalLikelihoods();
    double prob = 0.0;
    for (int c=0; c<numIllegalPats; c++)
        {
        double like = pi[space][0] * clP[0] + pi[space][1] * clP[1];
        prob += like;
        clP += 2;
        }

    return prob;
}

double Model::lnPriorRatio(int up, int down) {

    return theTree[up]->lnPrior() - theTree[down]->lnPrior();
}

void Model::mapCharacters(int space) {

    Tree* t = theTree[space];

    // pass down the tree initializing the conditional likelihoods
    for (int n=0; n<t->getNumberOfDownPassNodes(); n++)
        {
        Node* p = t->getDownPassNode(n);
        if (p->isLeaf() == false)
            {
            std::set<Node*> descendants = p->getDescendants();
            double* clP = p->exposeConditionalLikelihoods();
            std::vector<double*> clQ(descendants.size());
            std::vector<double**> ti(descendants.size());
            int k = 0;
            for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
                {
                clQ[k] = (*it)->exposeConditionalLikelihoods();
                ti[k] = (*it)->exposeTransitionProbability();
                k++;
                }
                
            for (int c=0; c<numGenes; c++)
                {
                for (int i=0; i<2; i++)
                    {
                    double product = 1.0;
                    for (int d=0; d<clQ.size(); d++)
                        {
                        double sum = 0.0;
                        for (int j=0; j<2; j++)
                            sum += ti[d][i][j] * clQ[d][j];
                        product *= sum;
                        }
                    clP[i] = product;
                    }
                clP += 2;
                for (int d=0; d<clQ.size(); d++)
                    clQ[d] += 2;
                }
                
            }
        }

    // pass up the tree, choosing a state at each interior node
    std::vector<int*> states(t->getNumNodes());
    for (int n=t->getNumberOfDownPassNodes()-1; n>=0; n--)
        {
        Node* p = t->getDownPassNode(n);
        states[p->getIndex()] = new int[numGenes];
        if (p == t->getRoot())
            {
            // ancestral node
            double* clP = p->exposeConditionalLikelihoods();
            for (int c=0; c<numGenes; c++)
                {
                double pr[2];
                pr[0] = clP[0] * pi[space][0];
                pr[1] = clP[1] * pi[space][1];
                double sum = pr[0] + pr[1];
                pr[0] /= sum;
                pr[1] /= sum;
                double u = ranPtr->uniformRv();
                if (u < pr[0])
                    states[p->getIndex()][c] = 0;
                else
                    states[p->getIndex()][c] = 1;
                clP += 2;
                }
            }
        else
            {
            // interior or tip node
            std::set<Node*> descendants = p->getDescendants();
            double* clP = p->exposeConditionalLikelihoods();
            double** tiP = p->exposeTransitionProbability();
            for (int c=0; c<numGenes; c++)
                {
                for (int i=0; i<2; i++)
                    {
                    clP[i] *= tiP[ states[p->getAncestor()->getIndex()][c] ][i];
                    }
                    
                double pr[2];
                pr[0] = clP[0] * pi[space][0];
                pr[1] = clP[1] * pi[space][1];
                double sum = pr[0] + pr[1];
                pr[0] /= sum;
                pr[1] /= sum;
                double u = ranPtr->uniformRv();
                if (u < pr[0])
                    states[p->getIndex()][c] = 0;
                else
                    states[p->getIndex()][c] = 1;
                    
                clP += 2;
                }
            }
        }
    
    // visit each branch and choose a history using rejection sampling
    double q[2][2];
    rateMatrix(space, q);
    for (int c=0; c<numGenes; c++)
        {
        CharacterChanges* changesSummary = characterChanges[c];
        int num01 = 0, num10 = 0;
        for (int n=0; n<t->getNumberOfDownPassNodes(); n++)
            {
            Node* p = t->getDownPassNode(n);
            if (p != t->getRoot())
                {
                BranchChanges* branchSummary = branchChanges[p->getIndex()];
                bool happyEnding = false;
                int n01 = 0, n10 = 0;
                do
                    {
                    double t = p->getBranchLength(), curT = 0.0;
                    int curState = states[ p->getAncestor()->getIndex() ][c];
                    int endState = states[ p->getIndex()                ][c];
                    double lambda = -q[curState][curState];
                    int numChanges = 0;
                    n01 = n10 = 0;
                    do
                        {
                        double u = ranPtr->uniformRv();
                        if (numChanges == 0 && curState != endState)
                            curT += -log(1.0 - u*(1.0-exp(-lambda*t))) / lambda;
                        else
                            curT += -log(u) / lambda;
                        if (curT < t)
                            {
                            if (curState == 0)
                                {
                                curState = 1;
                                n01++;
                                }
                            else
                                {
                                curState = 0;
                                n10++;
                                }
                            lambda = -q[curState][curState];
                            numChanges++;
                            }
                        } while (curT < t);
                    if (curState == endState)
                        happyEnding = true;
                    } while (happyEnding == false);
                branchSummary->updateChangesForSite(c, n01, n10, states[p->getIndex()][c]);
                num01 += n01;
                num10 += n10;
                }
            }
        changesSummary->updateChanges(num01, num10);
        }

    // free the changes vector
    for (int i=0; i<states.size(); i++)
        {
        if (states[i] != NULL)
            delete [] states[i];
        }

}

void Model::rateMatrix(int space, double q[2][2]) {

    q[0][1] = pi[space][1];
    q[1][0] = pi[space][0];
    q[0][0] = -q[0][1];
    q[1][1] = -q[1][1];
    double aveRate = pi[space][0] * q[0][1] + pi[space][1] * q[1][0];
    q[0][0] /= aveRate;
    q[0][1] /= aveRate;
    q[1][0] /= aveRate;
    q[1][1] /= aveRate;
}

bool Model::readTree(void) {

	// open the file
	std::ifstream treeStream(settingsPtr->getTreeFileName().c_str());
	if (!treeStream)
		{
		std::cerr << "Cannot open file \"" + settingsPtr->getTreeFileName() + "\"" << std::endl;
		return false;
		}
    
    std::string linestring = "";
	while( getline(treeStream, linestring).good() )
		{
		std::istringstream linestream(linestring);
        if (linestring != "")
            {
            theTree[0] = new Tree(linestring, settingsPtr, dataPtr, this, ranPtr, branchLengthPriorParm, 0);
            }
        }

	treeStream.close();
    theTree[0]->getDownPassSequence();
    
    // mark clades
    for (int cld=0; cld<dataPtr->getNumClades(); cld++)
        {
        Clades* theClade = dataPtr->getClade(cld);
        theTree[0]->setAllFlagsTo(false);
        
        Node* nde = theTree[0]->getNodeWithTaxon(theClade->getTaxon(0));
        if (nde == NULL)
            Msg::error("Problem finding first taxon in clade");
        nde->setFlag(true);

        nde = theTree[0]->getNodeWithTaxon(theClade->getTaxon(1));
        if (nde == NULL)
            Msg::error("Problem finding second taxon in clade");
        nde->setFlag(true);
        
        for (int n=0; n<theTree[0]->getNumberOfDownPassNodes(); n++)
            {
            nde = theTree[0]->getDownPassNode(n);
            int npft = nde->numberOfDescendantsWithFlag(true);
            if (npft == 1)
                nde->setFlag(true);
            else if (npft >= 2)
                {
                nde->setFlag(true);
                nde->setName(theClade->getName());
                nodesOfInterest.push_back(nde);
                break;
                }
            }
        }
    
    return true;
}

void Model::reject(void) {

    *theTree[1] = *theTree[0];
    pi[1][0] = pi[0][0];
    pi[1][1] = pi[0][1];
}

double Model::update(int space) {

    double u = ranPtr->uniformRv();
    double sum = 0.0;
    int whichMove = 0;
    for (int i=0; i<proposalProbs.size(); i++)
        {
        sum += proposalProbs[i];
        if (u < sum)
            {
            whichMove = i;
            break;
            }
        }
    if (whichMove == 0)
        return theTree[space]->update();
    else
        return updatePi(space);
}

double Model::updatePi(int space) {

    double alpha0 = 200.0;
    std::vector<double> aForward(2);
    std::vector<double> aReverse(2);
    std::vector<double> oldFreqs(2);
    std::vector<double> newFreqs(2);
    for (int i=0; i<2; i++)
        {
	oldFreqs[i] = pi[space][i];
	aForward[i] = pi[space][i] * alpha0;
	}
    ranPtr->dirichletRv(aForward, newFreqs);
    for (int i=0; i<2; i++)
        {
	aReverse[i] = newFreqs[i] * alpha0;
	pi[space][i] = newFreqs[i];
        }
    
    //std::cout << "Model::updatePi(" << space << ")\tpi[0]\t" << pi[space][0] << "\tpi[1]\t" << pi[space][1] << std::endl;

    std::vector<Node*> theTreeNodes = theTree[space]->exposeNodes();
    for (int i=0; i<theTreeNodes.size(); i++)
        {
        if ( theTreeNodes[i] != theTree[space]->getRoot() )
            {
            double x = theTreeNodes[i]->getBranchProportion();
            //std::cout << "Model::updatePi calls setBranchProportion(" << x << ")";
	    theTreeNodes[i]->setBranchProportion(x);
            }
        }
    
    return ranPtr->lnDirichletPdf(aReverse, oldFreqs) - ranPtr->lnDirichletPdf(aForward, newFreqs);
}

double* Model::getPi(int space) {

    //std::cout << "Model::getPi\tspace\t" << space << "\tpi[0]\t" << pi[space][0] << "\tpi[1]\t" << pi[space][1] << std::endl;
    return pi[space];
}
