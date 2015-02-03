#ifndef Tree_H
#define Tree_H

#include <sstream>
#include <string>
#include <vector>

class Data;
class MbRandom;
class Model;
class Node;
class Settings;

class Tree {
    
	public:
                            Tree(std::string treeStr, Settings* sp, Data* dp, Model* mp, MbRandom* rp, double lam, int space);
                            Tree(Tree& t, int space);
                           ~Tree(void);
        Tree&               operator=(Tree& t);
        std::vector<Node*>  exposeNodes(void) { return nodes; }
        Node*               getDownPassNode(int idx) { return downPassSequence[idx]; }
        int                 getNumberOfDownPassNodes(void) { return (int)downPassSequence.size(); }
        void                getDownPassSequence(void);
        int                 getNumNodes(void) { return (int)nodes.size(); }
        Node*               getRoot(void) { return root; }
        std::string         getTreeString(void);
        double              getTreeLength(void) { return treeLength; }
        void                markPathToRootFromNode(Node* p);
        double              lnPrior(void);
        void                print(void);
        void                setAllFlagsTo(bool tf);
        void                showNodeInfo(void);
        Node*               getNodeWithTaxon(std::string s);
        Node*               getNodeWithIndex(int idx) { return nodes[idx]; }
        double              update(void);
        void                setMySpace(int x) { mySpace = x; }
        int                 getMySpace(void) { return mySpace; }
        void                setTreeLength(double x) { treeLength = x; }

    private:
        void                clone(Tree& t);
        void                deleteNodes(void);
        void                passDn(Node* p);
        void                initializePartitions(void);
        double              updateBranchProportions(void);
        double              updateTreeLength(void);
        Data*               dataPtr;
        MbRandom*           ranPtr;
        Model*              modelPtr;
        Settings*           settingsPtr;
        void                interpretTreeString(std::string treeStr);
        void                showNodes(Node* p, int indent);
        void                writeTree(Node* p, std::stringstream& ss);
        std::vector<Node*>  nodes;
        std::vector<Node*>  downPassSequence;
        Node*               root;
        int                 mySpace;
        double              treeLength;
        double              lambda;
};


#endif
