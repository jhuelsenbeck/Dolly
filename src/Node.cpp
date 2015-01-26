#include "Model.h"
#include "Msg.h"
#include "Node.h"
#include "Tree.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>



Node::Node(Model* mp, Tree* t, int nc) {

    modelPtr = mp;
    myTree   = t;
    index    = 0;
    ancestor = NULL;
    name     = "";
    leafNode = false;
    numChar  = nc;
    if (numChar <= 0)
        Msg::error("Too few characters");
    
    transitionProbabilities = new double*[2];
    transitionProbabilities[0] = new double[4];
    for (int i=1; i<2; i++)
        transitionProbabilities[i] = transitionProbabilities[i-1] + 2;
    for (int i=0; i<2; i++)
        for (int j=0; j<2; j++)
            transitionProbabilities[i][j] = 0.0;
    
    cls = new double[2*numChar];
    for (int i=0; i<2*numChar; i++)
        cls[i] = 0.0;

    setBranchProportion(1.0);
    descendants.clear();
}

Node::Node(void) {

}

Node::~Node(void) {

    delete [] transitionProbabilities[0];
    delete [] transitionProbabilities;
    delete [] cls;
}

double Node::getBranchLength(void) {

    return branchProportion * myTree->getTreeLength();
}

std::string Node::getDescendantsString(void) {

    std::string str = "";
    char tempCharStr[20];
    for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
        if (it == descendants.begin())
            sprintf(tempCharStr, "%d", (*it)->getIndex());
        else
            sprintf(tempCharStr, ",%d", (*it)->getIndex());
        str += tempCharStr;
        }
    
    return str;
}

bool Node::isParentFlag(bool tf) {

    for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
        if ( (*it)->getFlag() == tf )
            return true;
        }
    return false;
}

int Node::numberOfDescendantsWithFlag(bool tf) {

    int x = 0;
    for (std::set<Node*>::iterator it = descendants.begin(); it != descendants.end(); it++)
        {
        if ( (*it)->getFlag() == tf )
            x++;
        }
    return x;
}

void Node::print(void) {

    std::cout << "address          = " << this << std::endl;
    std::cout << "index            = " << index << std::endl;
    std::cout << "memoryIdx        = " << memoryIdx << std::endl;
    std::cout << "branchProportion = " << branchProportion << std::endl;
    std::cout << "name             = " << name << std::endl;
    std::cout << "modelPtr         = " << modelPtr << std::endl;
    std::cout << "leafNode         = " << leafNode << std::endl;
    std::cout << "flag             = " << flag << std::endl;
    if (ancestor != NULL)
        std::cout << "ancestor     = " << ancestor->index << std::endl;
    else
        std::cout << "ancestor     = " << "NULL" << std::endl;
    std::cout << "descendants  = ";
    for (std::set<Node*>::iterator it=descendants.begin(); it != descendants.end(); it++)
        std::cout << (*it)->index << " ";
    std::cout << std::endl;
    std::cout << "myPartition  = ";
    for (std::set<std::string>::iterator it=myPartition.begin(); it != myPartition.end(); it++)
        std::cout << *it << " ";
    std::cout << std::endl;
    std::cout << "               " << transitionProbabilities[0][0] << " " << transitionProbabilities[0][1] << std::endl;
    std::cout << "               " << transitionProbabilities[1][0] << " " << transitionProbabilities[1][1] << std::endl;
}

void Node::printConditionalLikelihoods(int prec) {

    std::cout << std::fixed << std::setprecision(prec);
    double* p = cls;
    for (int i=0; i<numChar; i++)
        {
        std::cout << p[0] << " " << p[1] << std::endl;
        p += 2;
        }
    std::cout << transitionProbabilities[0][0] << " " << transitionProbabilities[0][1] << std::endl;
    std::cout << transitionProbabilities[1][0] << " " << transitionProbabilities[1][1] << std::endl;
}

void Node::setBranchProportion(double x) {

    branchProportion = x;
    double branchLength = getBranchLength();
    double* stationaryFrequencies = modelPtr->getPi(myTree->getMySpace());
    double u = 1.0 / (2.0 * stationaryFrequencies[0] * stationaryFrequencies[1]);
    double expPart = exp(-u*branchLength);
    transitionProbabilities[0][0] = stationaryFrequencies[0] + stationaryFrequencies[1]*expPart;
    transitionProbabilities[0][1] = stationaryFrequencies[1] - stationaryFrequencies[1]*expPart;
    transitionProbabilities[1][0] = stationaryFrequencies[0] - stationaryFrequencies[0]*expPart;
    transitionProbabilities[1][1] = stationaryFrequencies[1] + stationaryFrequencies[0]*expPart;
#   if 0
    std::cout << std::fixed << std::setprecision(4);
    for (int i=0; i<2; i++)
        {
        for (int j=0; j<2; j++)
            {
            std::cout << transitionProbabilities[i][j] << " ";
            }
        std::cout << std::endl;
        }
#   endif
}
