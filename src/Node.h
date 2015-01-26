#ifndef Node_H
#define Node_H

#include <set>
#include <string>
class Model;
class Tree;

class Node {

    friend class                Tree;
    
	public:
                                Node(Model* mp, Tree* t, int nc);
                               ~Node(void);
        void                    addDescendant(Node* p) { descendants.insert(p); }
        void                    addTaxonToPartition(std::string s) { myPartition.insert(s); }
        double*                 exposeConditionalLikelihoods(void) { return cls; }
        double*                 exposeDummyConditionalLikelihoods(void) { return dummyCls; }
        double**                exposeTransitionProbability(void) { return transitionProbabilities; }
        Node*                   getAncestor(void) { return ancestor; }
        double                  getBranchLength(void);
        double                  getBranchProportion(void) { return branchProportion; }
        std::set<Node*>&        getDescendants(void) { return descendants; }
        std::set<std::string>   getPartition(void) { return myPartition; }
        std::string             getDescendantsString(void);
        bool                    getFlag(void) { return flag; }
        int                     getIndex(void) { return index; }
        size_t                  getMemoryIdx(void) { return memoryIdx; }
        std::string             getName(void) { return name; }
        int                     getNumberOfDescendants(void) { return (int)descendants.size(); }
        bool                    isLeaf(void) { return leafNode; }
        bool                    isParentFlag(bool tf);
        int                     numberOfDescendantsWithFlag(bool tf);
        void                    print(void);
        void                    printConditionalLikelihoods(int prec);
        void                    removeDescendant(Node* p) { descendants.erase(p); }
        void                    removeDescendants(void) { descendants.clear(); }
        void                    setBranchProportion(double x);
        void                    setFlag(bool tf) { flag = tf; }
        void                    setIndex(int x ) { index = x; }
        void                    setLeaf(bool tf) { leafNode = tf; }
        void                    setName(std::string s) { name = s; }
        void                    setMemoryIdx(size_t x) { memoryIdx = x; }
        void                    setAncestor(Node* p) { ancestor = p; }

    private:
                                Node(void);
        Model*                  modelPtr;
        Tree*                   myTree;
        int                     index;
        size_t                  memoryIdx;
        double                  branchProportion;
        std::string             name;
        bool                    leafNode;
        Node*                   ancestor;
        std::set<Node*>         descendants;
        double**                transitionProbabilities;
        std::set<std::string>   myPartition;
        bool                    flag;
        int                     numChar;
        double*                 cls;
        double                  dummyCls[4];
};


#endif
