#include <iostream>
#include "Data.h"
#include "MbRandom.h"
#include "Mcmc.h"
#include "Model.h"
#include "Settings.h"
// #include "Util.h"

void printHeader(void);


int main(int argc, char* argv[]) {

    // print a header
    printHeader();

    // instantiate the random number object
    MbRandom myRandom;
    
    // read user settings
    Settings mySettings( argc, argv );
    
    // read the data
    Data myData( mySettings.getInputFileName(), mySettings.getGenePresenceProbabilitiesFileName(), mySettings.getUseGenePresenceProbs() );
    
    // set up the phylogenetic model
    Model myModel( &mySettings, &myData, &myRandom );
    
    // run a Markov chain to integrate over branch-length and model parameter combinations
    Mcmc myMcmc( &mySettings, &myModel, &myRandom );
    
    // calculate the probabilities at each interior node on the tree
    //myModel.calculateStateProbs(0);


    return 0;
}

void printHeader(void) {

    std::cout << std::endl;
	std::cout << "   Dolly" << std::endl;
	std::cout << "   John P. Huelsenbeck & Daniel J. Richter" << std::endl;
	std::cout << "   University of California, Berkeley" << std::endl;
    std::cout << std::endl;
}
