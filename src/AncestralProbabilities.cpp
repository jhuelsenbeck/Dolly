#include "AncestralProbabilities.h"



AncestralProbabilities::AncestralProbabilities(Node* np, int nc) {

    numChar = nc;
    myNode  = np;
    probabilityOfOne = new double[numChar];
    for (int i=0; i<numChar; i++)
        probabilityOfOne[i] = 0.0;

}

AncestralProbabilities::~AncestralProbabilities(void) {

    delete [] probabilityOfOne;
}
