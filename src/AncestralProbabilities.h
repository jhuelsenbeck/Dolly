#ifndef AncestralProbabilities_H
#define AncestralProbabilities_H

class Node;

class AncestralProbabilities {

	public:
                                AncestralProbabilities(Node* np, int nc);
                               ~AncestralProbabilities(void);

    private:
                                AncestralProbabilities(void) { };
        double*                 probabilityOfOne;
        int                     numChar;
        Node*                   myNode;
};


#endif
