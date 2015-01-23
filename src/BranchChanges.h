#ifndef BranchChanges_H
#define BranchChanges_H


#include <string>
class Node;

class BranchChanges {

	public:
                                    BranchChanges(int nc, int idx, std::string nme);
                                   ~BranchChanges(void);
        std::string                 getName(void) { return name; }
        double                      getAverageNumGains(int idx);
        double                      getAverageNumLosses(int idx);
        double                      getProbabilityOne(int idx);
        void                        updateChangesForSite(int site, int n01, int n10, int endState);

    private:
        std::string                 name;
        int                         myNodeIndex;
        int                         numCharacters;
        int*                        numSamples;
        int*                        numGains;
        int*                        numLosses;
        int*                        numOneAtNode;
};


#endif
