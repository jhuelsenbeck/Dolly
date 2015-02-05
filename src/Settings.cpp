#include "Msg.h"
#include "Settings.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>


Settings::Settings(int argc, char *argv[]) {

#   if 0
	// set up fake command-line argument string
	char* cmdString[24];
	cmdString[ 0] = (char *)"dollo";
	cmdString[ 1] = (char *)"-input_file";
	cmdString[ 2] = (char *)"/Users/johnh/Dolly/example_data/single_site/gene_presence_absence";
  //cmdString[ 2] = (char *)"/Users/johnh/Desktop/DolloPlusData/states";
	cmdString[ 3] = (char *)"-tree_file";
	cmdString[ 4] = (char *)"/Users/johnh/Dolly/example_data/single_site/with_branch_lengths.tree";
  //cmdString[ 4] = (char *)"/Users/johnh/Desktop/DolloPlusData/unrooted_tree_no_branch_lengths";
	cmdString[ 5] = (char *)"-error_file";
	cmdString[ 6] = (char *)"/Users/johnh/Desktop/DolloPlusData/probabilities";
	cmdString[ 7] = (char *)"-output_file";
	cmdString[ 8] = (char *)"/Users/johnh/Desktop/DolloPlusResults/out";
	cmdString[ 9] = (char *)"-length";
	cmdString[10] = (char *)"1000";
	cmdString[11] = (char *)"-burn";
	cmdString[12] = (char *)"0";
	cmdString[13] = (char *)"-print_freq";
	cmdString[14] = (char *)"1";
	cmdString[15] = (char *)"-sample_freq";
	cmdString[16] = (char *)"10";
	cmdString[17] = (char *)"-use_errors";
	cmdString[18] = (char *)"no";
	cmdString[19] = (char *)"-brlen_prior";
	cmdString[20] = (char *)"40.0";
	cmdString[21] = (char *)"-user_brlens";
	cmdString[22] = (char *)"yes";
	argc = 23;
	argv = cmdString;
#   endif

    // default values
    chainLength                    = 1000;
    burnIn                         = 0;
    useErrorProbs                  = false;
    inputFileName                  = "";
    errorProbabilitiesFileName     = "";
    treeFileName                   = "";
    outputFileName                 = "";
    printFrequency                 = 10;
    sampleFrequency                = 10;
    branchLengthPrior              = 40.0;
    fixBranchProportionsToUserTree = false;
    
	if (argc > 1)
		{
		if (argc % 2 == 0)
			printUsage();
        
		/* read the command-line arguments */
		std::string status = "none";
		for (int i=1; i<argc; i++)
			{
			std::string cmd = argv[i];
            //std::cout << cmd << std::endl;
			if (status == "none")
				{
				/* read the parameter specifier */
				if ( cmd == "-input_file" )
					status = "input_file";
				else if ( cmd == "-output_file" )
					status = "output_file";
				else if ( cmd == "-tree_file" )
					status = "tree_file";
				else if ( cmd == "-error_file" )
					status = "error_file";
				else if ( cmd == "-use_errors" )
					status = "use_errors";
				else if ( cmd == "-length" )
					status = "length";
				else if ( cmd == "-burn" )
					status = "burn";
				else if ( cmd == "-print_freq" )
					status = "print_freq";
				else if ( cmd == "-sample_freq" )
					status = "sample_freq";
				else if ( cmd == "-brlen_prior" )
					status = "brlen_prior";
				else if ( cmd == "-user_brlens" )
					status = "user_brlens";
				else
					{
					std::cerr << "Could not interpret option \"" << cmd << "\"." << std::endl;
					exit(1);
					}
				}
			else
				{
				/* read the parameter */
				if ( status == "input_file" )
					{
					inputFileName = argv[i];
					}
				else if ( status == "tree_file" )
					{
					treeFileName = argv[i];
					}
				else if ( status == "error_file" )
					{
					errorProbabilitiesFileName = argv[i];
					}
				else if ( status == "use_errors" )
					{
                    if ( (argv[i][0] == 'Y' || argv[i][0] == 'y') &&
                         (argv[i][1] == 'E' || argv[i][1] == 'e') &&
                         (argv[i][2] == 'S' || argv[i][2] == 's') )
                        useErrorProbs = true;
                    else if ( (argv[i][0] == 'N' || argv[i][0] == 'n') &&
                              (argv[i][1] == 'O' || argv[i][1] == 'o') )
                        useErrorProbs = false;
                    else
                        Msg::error("Unknown option for use_errors");
					}
				else if ( status == "output_file" )
					{
					outputFileName = argv[i];
					}
				else if ( status == "length" )
					{
					chainLength = atoi(argv[i]);
					}
				else if ( status == "burn" )
					{
					burnIn = atoi(argv[i]);
					}
				else if ( status == "print_freq" )
					{
					printFrequency = atoi(argv[i]);
					}
				else if ( status == "sample_freq" )
					{
					sampleFrequency = atoi(argv[i]);
					}
				else if ( status == "brlen_prior" )
					{
					branchLengthPrior = (double)atof(argv[i]);
					}
				else if ( status == "user_brlens" )
					{
                    std::string str = argv[i];
                    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
                    if (str[0] == 'y' && str[1] == 'e' && str[2] == 's')
                        fixBranchProportionsToUserTree = true;
                    else if (str[0] == 'n' && str[1] == 'o')
                        fixBranchProportionsToUserTree = false;
                    else
                        Msg::error("Unknown option for command \"user_brlens\"");
					}
				else
					{
					Msg::error("Unknown status reading command line information");
					}
				status = "none";
				}
			}
		}
	else
		{
		printUsage();
		}

    print();
}

void Settings::print(void) {
    
	std::cout << "   Analysis information:" << std::endl;
	std::cout << "   * Input file name                     = " << inputFileName                  << std::endl;
	std::cout << "   * Tree file name                      = " << treeFileName                   << std::endl;
	std::cout << "   * Output file name                    = " << outputFileName                 << std::endl;
	std::cout << "   * Error file name                     = " << errorProbabilitiesFileName     << std::endl;
	std::cout << "   * Use error probabilities             = " << (useErrorProbs == 1 ? "yes" : "no") << std::endl;
	std::cout << "   * Chain length                        = " << chainLength                    << std::endl;
	std::cout << "   * Burn In                             = " << burnIn                         << std::endl;
	std::cout << "   * Print frequency                     = " << printFrequency                 << std::endl;
	std::cout << "   * Parameter sample frequency          = " << sampleFrequency                << std::endl;
	std::cout << "   * Branch length prior parm            = " << branchLengthPrior              << std::endl;
	std::cout << "   * Fix branch length prop to user tree = " << (fixBranchProportionsToUserTree == 1 ? "yes" : "no") << std::endl;
	std::cout << std::endl;
}

void Settings::printUsage(void) {
    
	std::cout << "Usage:" << std::endl;
	std::cout << "   -input_file <FILE NAME>  : Gene presence/absence file name [required]" << std::endl;
	std::cout << "   -tree_file <FILE NAME>   : Tree file name in Newick format [required]" << std::endl;
	std::cout << "   -output_file <FILE NAME> : Output file name [required]" << std::endl;
	std::cout << "   -error_file <FILE NAME>  : File containing probabilities for gene presence [optional]" << std::endl;
        std::cout << "   -use_errors <yes/no>     : Whether to use probabilities for gene presence [default " << (useErrorProbs == 1 ? "yes" : "no") << "]" << std::endl;
	std::cout << "   -length <NUMBER>         : Number of MCMC cycles [default " << chainLength << "]" << std::endl;
	std::cout << "   -burn <NUMBER>           : Number of MCMC cycles to discard (burn) [default " << burnIn << "]" << std::endl;
	std::cout << "   -print_freq <NUMBER>     : Frequency with which information is printed to the screen [default " << printFrequency << "]" << std::endl;
	std::cout << "   -sample_freq <NUMBER>    : Frequency with which information is printed to files [default " << sampleFrequency << "]" << std::endl;
	std::cout << "   -brlen_prior <NUMBER>    : Branch length prior parameter [default " << branchLengthPrior << "]" << std::endl;
	std::cout << "   -user_brlens <yes/no>    : Fix branch length proportions to user tree [default " << (fixBranchProportionsToUserTree == 1 ? "yes" : "no") << "]" << std::endl;
    
	std::cout << std::endl;
	std::cout << "Example:" << std::endl;
	std::cout << "   ./dolly -input_file <input file> -tree_file <tree file> -output_file <output file>" << std::endl;
	exit(1);
}
