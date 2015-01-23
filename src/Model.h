#ifndef Model_H
#define Model_H

#include <vector>
class AncestralProbabilities;
class BranchChanges;
class CharacterChanges;
class Data;
class MbRandom;
class Node;
class Settings;
class Tree;

class Model {

	public:
                                                Model(Settings* sp, Data* dp, MbRandom* rp);
                                               ~Model(void);
        void                                    accept(void);
        double                                  averageBranchLength(int space);
        void                                    calculateStateProbs(int space);
        BranchChanges*                          getBranchChangesForNode(int nde) { return branchChanges[nde]; }
        double*                                 getPi(int space) { return pi[space]; }
        Tree*                                   getTree(int space) { return theTree[space]; }
        int                                     getNumGenes(void) { return numGenes; }
        double                                  lnLikelihood(void);
        double                                  lnLikelihood(int space);
        double                                  lnPriorRatio(int up, int down);
        void                                    mapCharacters(int space);
        void                                    reject(void);
        double                                  update(int space);
        CharacterChanges*                       getCharacterChange(int idx) { return characterChanges[idx]; }

    private:
        void                                    initializeConditionalLikelihood(void);
        void                                    calcConditionalLikelihoodsDown(int space);
        void                                    calculateLikelihoodsDownTree(int space);
        void                                    calcConditionalLikelihoodsUpToNode(int space, Node* theNode);
        void                                    calculateStateProbabilitiesAtNode(int space, Node* theNode);
        bool                                    readTree(void);
        double                                  updatePi(int space);
        double                                  probInvariant(int space);
        void                                    rateMatrix(int space, double q[2][2]);
        Data*                                   dataPtr;
        Settings*                               settingsPtr;
        MbRandom*                               ranPtr;
        Tree*                                   theTree[2];
        double                                  pi[2][2];
        std::vector<double>                     alpha;
        std::vector<AncestralProbabilities*>    ancestralProbabilities;
        int                                     numGenes;
        int                                     numNodes;
        int                                     numStates;
        std::vector<Node*>                      nodesOfInterest;
        std::vector<double>                     proposalRates;
        std::vector<double>                     proposalProbs;
        std::vector<CharacterChanges*>          characterChanges;
        std::vector<BranchChanges*>             branchChanges;
};


#endif
