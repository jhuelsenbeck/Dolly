#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include "Data.h"
#include "MbRandom.h"
#include "Msg.h"
#include "Node.h"
#include "Settings.h"
#include "Tree.h"


Tree::Tree(std::string treeStr, Settings *sp, Data* dp, Model* mp, MbRandom* rp, double lam, int space) {

    dataPtr     = dp;
    modelPtr    = mp;
    ranPtr      = rp;
    settingsPtr = sp;
    mySpace     = space;
    lambda      = lam;
    interpretTreeString(treeStr);
    //print();
}

Tree::Tree(Tree& t, int space) {

    mySpace = space;
    if ( this != &t)
        clone(t);
}

Tree::~Tree(void) {

    deleteNodes();
}

Tree& Tree::operator=(Tree& t) {

    if ( this != &t)
        clone(t);
    return *this;
}

void Tree::clone(Tree& t) {

    dataPtr     = t.dataPtr;
    ranPtr      = t.ranPtr;
    modelPtr    = t.modelPtr;
    settingsPtr = t.settingsPtr;
    treeLength  = t.treeLength;
    lambda      = t.lambda;
    
    if (nodes.size() != t.nodes.size())
        {
        deleteNodes();
        for (int i=0; i<t.nodes.size(); i++)
            {
            Node* newNode = new Node(modelPtr, this, dataPtr->getNumChar());
            newNode->setMemoryIdx(i);
            nodes.push_back( newNode );
            }
        }
    
    for (int n=0; n<nodes.size(); n++)
        {
        Node* sourceNode = t.nodes[n];
        Node* myNode     = nodes[n];
        
        if (sourceNode == t.root)
            root = myNode;

        myNode->modelPtr         = sourceNode->modelPtr;
        myNode->index            = sourceNode->index;
        myNode->branchProportion = sourceNode->branchProportion;
        myNode->name             = sourceNode->name;
        myNode->leafNode         = sourceNode->leafNode;
        myNode->myPartition      = sourceNode->myPartition;
        myNode->flag             = sourceNode->flag;
        for (int i=0; i<2; i++)
            for (int j=0; j<2; j++)
                myNode->transitionProbabilities[i][j] = sourceNode->transitionProbabilities[i][j];
        if (sourceNode->ancestor != NULL)
            myNode->ancestor = nodes[ sourceNode->ancestor->memoryIdx ];
        myNode->removeDescendants();
        std::set<Node*> descs = sourceNode->getDescendants();
        for (std::set<Node*>::iterator it = descs.begin(); it != descs.end(); it++)
            {
            Node* nde = nodes[ (*it)->memoryIdx ];
            myNode->addDescendant(nde);
            }
        }
    
    getDownPassSequence();
}

void Tree::deleteNodes(void) {

    for (std::vector<Node*>::iterator it=nodes.begin(); it != nodes.end(); it++)
        delete (*it);
    nodes.clear();
}

Node* Tree::getNodeWithTaxon(std::string s) {

    for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        {
        if ( (*it)->getName() == s )
            return (*it);
        }
    return NULL;
}

void Tree::getDownPassSequence(void) {

    downPassSequence.clear();
    passDn(root);
}

void Tree::initializePartitions(void) {

    for (int n=0; n<getNumberOfDownPassNodes(); n++)
        {
        Node* p = getDownPassNode(n);
        if (p->isLeaf() == true)
            {
            p->addTaxonToPartition(p->getName());
            }
        else
            {
            std::set<Node*> des = p->getDescendants();
            for (std::set<Node*>::iterator it = des.begin(); it != des.end(); it++)
                {
                std::set<std::string> part = (*it)->getPartition();
                for (std::set<std::string>::iterator it2 = part.begin(); it2 != part.end(); it2++)
                    {
                    p->addTaxonToPartition(*it2);
                    }
                }
            }
        }
}

void Tree::interpretTreeString(std::string treeStr) {

    std::vector<Node*> interiorNodes;
    std::vector<Node*> tipNodes(dataPtr->getNumTaxa());
    
    std::string nameStr = "";
    Node* p = NULL;
    int interiorNodeIdx = dataPtr->getNumTaxa();
    bool readingTaxon = true;
    bool foundUserBrlens = false;
    for (int i=0; i<treeStr.size(); i++)
        {
        char c = treeStr[i];
        if (c == '(')
            {
            Node* newNode = new Node(modelPtr, this, dataPtr->getNumChar());
            newNode->setIndex(interiorNodeIdx++);
            if (p == NULL)
                {
                root = newNode;
                }
            else
                {
                p->addDescendant(newNode);
                newNode->setAncestor(p);
                }
            interiorNodes.push_back(newNode);
            p = newNode;
            }
        else if (c == ')')
            {
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            }
        else if (c == ',')
            {
            if (p->getAncestor() != NULL)
                p = p->getAncestor();
            }
        else if (c == ';')
            {
            }
        else if (c == ':')
            {
            readingTaxon = false;
            }
        else if (c != ' ')
            {
            nameStr += c;
            if ( treeStr[i+1] == ')' || treeStr[i+1] == '(' ||
                 treeStr[i+1] == ',' || treeStr[i+1] == ';' || treeStr[i+1] == ':' )
                {
                if (readingTaxon == true)
                    {
                    int taxonIdx = dataPtr->getTaxonIndex(nameStr);
                    if (taxonIdx == -1)
                        Msg::error("Problem finding taxon " + nameStr);
                    Node* newNode = new Node(modelPtr, this, dataPtr->getNumChar());
                    newNode->setLeaf(true);
                    newNode->setName(nameStr);
                    newNode->setIndex(taxonIdx);
                    if (p == NULL)
                        {
                        root = newNode;
                        }
                    else
                        {
                        p->addDescendant(newNode);
                        newNode->setAncestor(p);
                        }
                    tipNodes[taxonIdx] = newNode;
                    p = newNode;
                    }
                else
                    {
                    double x;
                    std::istringstream buf(nameStr);
                    buf >> x;
                    p->setBranchProportion(x);
                    readingTaxon = true;
                    foundUserBrlens = true;
                    }
                nameStr = "";
                }
            }
        }
    
    // add the tip and interior nodes to the class instance variable, nodes
    for (int i=0; i<tipNodes.size(); i++)
        nodes.push_back(tipNodes[i]);
    for (int i=0; i<interiorNodes.size(); i++)
        nodes.push_back(interiorNodes[i]);
        
    // set memory indices
    for (int i=0; i<nodes.size(); i++)
        nodes[i]->setMemoryIdx(i);
    
    //
    if (settingsPtr->getFixBranchProportionsToUserTree() == false || foundUserBrlens == false)
        {
        for (int i=0; i<nodes.size(); i++)
            {
            if (nodes[i] != root)
                nodes[i]->setBranchProportion(1.0);
            }
        }
    
    // renormalize branch proportions to sum to one
    double sum = 0.0;
    for (int i=0; i<nodes.size(); i++)
        {
        if (nodes[i] != root)
            sum += nodes[i]->getBranchProportion();
        }
    treeLength = 0.0;
    std::cout << "lambda=" << lambda << " tl=" << treeLength << std::endl;
    for (int i=0; i<nodes.size(); i++)
        {
        if (nodes[i] != root)
            {
            double x = nodes[i]->getBranchProportion();
            nodes[i]->setBranchProportion(x / sum);
            treeLength += ranPtr->exponentialRv(lambda);
            }
        }
}

void Tree::markPathToRootFromNode(Node* p) {

    setAllFlagsTo(false);
    Node* q = p;
    q->setFlag(true);
    while (q->getAncestor() != NULL)
        {
        q->setFlag(true);
        q = q->getAncestor();
        }
    q->setFlag(true);
}

double Tree::lnPrior(void) {

    double lnP = 0.0;
    for (int i=0; i<nodes.size(); i++)
        {
        Node* p = nodes[i];
        if (p->ancestor != NULL)
            {
            lnP += (log(lambda) - lambda * p->getBranchLength());
            //std::cout << p->index << " v=" << p->getBranchLength() << std::endl;
            }
        }
    return lnP;
}

void Tree::passDn(Node* p) {

    if (p != NULL)
        {
        std::set<Node*> descendants = p->getDescendants();
        for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
            passDn(*it);
        downPassSequence.push_back(p);
        }
}

void Tree::print(void) {

    showNodes(root, 2);
}

void Tree::setAllFlagsTo(bool tf) {

    for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
        (*it)->setFlag(tf);
}

void Tree::showNodes(Node* p, int indent) {

    if (p != NULL)
        {
        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex() << " (" << p->getDescendantsString() << ") ";
        std::cout << p->getName();
        if (p == root)
            std::cout << " <- Root";
        std::cout << std::endl;
        std::set<Node*> descendantSet = p->getDescendants();
        for (std::set<Node*>::iterator it = descendantSet.begin(); it != descendantSet.end(); it++)
            showNodes(*it, indent+2);
        }
}

void Tree::showNodeInfo(void) {

    for (int i=0; i<nodes.size(); i++)
        nodes[i]->print();
}

double Tree::update(void) {

    double updateBranchProportionsWithProb = 0.8;
    if (settingsPtr->getFixBranchProportionsToUserTree() == true)
        updateBranchProportionsWithProb = 0.0;

    double u = ranPtr->uniformRv();
    if (u < updateBranchProportionsWithProb)
        return updateBranchProportions();
    else
        return updateTreeLength();
}

double Tree::updateBranchProportions(void) {

    // pick a node to update
    Node* nde = NULL;
    do
        {
        nde = nodes[(int)(nodes.size()*ranPtr->uniformRv())];
        } while(nde == root);

    // initialize some variables
    double alpha0 = 1000.0;
	std::vector<double> aForward(2);
	std::vector<double> aReverse(2);
	std::vector<double> oldProps(2);
    std::vector<double> newProps(2);

    // propose new values
    oldProps[0] = nde->getBranchProportion();
    oldProps[1] = 1.0 - oldProps[0];
    aForward[0] = oldProps[0] * alpha0;
    aForward[1] = oldProps[1] * alpha0;
	ranPtr->dirichletRv(aForward, newProps);
    aReverse[0] = newProps[0] * alpha0;
    aReverse[1] = newProps[1] * alpha0;
    
    // set the branch proportions (and update the branch lengths)
    double f = newProps[1]/oldProps[1];
    int n = 0;
    for (int i=0; i<nodes.size(); i++)
        {
        Node* tn = nodes[i];
        if (tn != root)
            {
            if ( tn == nde )
                tn->setBranchProportion(newProps[0]);
            else
                {
                double x = tn->getBranchProportion();
                tn->setBranchProportion(x * f);
                n++;
                }
            }
        }
    
    double lnProposalRatio = ranPtr->lnDirichletPdf(aReverse, oldProps) - ranPtr->lnDirichletPdf(aForward, newProps);
    lnProposalRatio += n * log(f); // Jacobian
    return lnProposalRatio;
}

double Tree::updateTreeLength(void) {

    double tuning = log(1.1);
    double oldL = treeLength;
    double newL = oldL * exp( tuning*(ranPtr->uniformRv()-0.5) );
    setTreeLength(newL);

    for (int i=0; i<nodes.size(); i++)
        {
        Node* tn = nodes[i];
        if (tn != root)
            {
            double x = tn->getBranchProportion();
            tn->setBranchProportion(x);
            }
        }
    
    return (log(newL) - log(oldL));
}

std::string Tree::getTreeString(void) {

    std::stringstream newickStream;
    writeTree(root, newickStream);
    std::string newickString = newickStream.str();
    return newickString;
}

void Tree::writeTree(Node* p, std::stringstream& ss) {

	if (p != NULL)
		{
		if (p->isLeaf() == true)
			{
			ss << p->getName() << ":" << std::fixed << std::setprecision(6) << p->getBranchLength();
			}
		else
			{
            ss << "(";
            std::set<Node*> myDescendants = p->getDescendants();
            int i = 0;
            for (std::set<Node*>::iterator it = myDescendants.begin(); it != myDescendants.end(); it++)
                {
                writeTree(*it, ss);
                if (i + 1 != myDescendants.size())
                    ss << ",";
                i++;
                }
            if (p == root)
                ss << ");";
            else
                ss << "):" << std::fixed << std::setprecision(6) << p->getBranchLength();
            }
		}
}


