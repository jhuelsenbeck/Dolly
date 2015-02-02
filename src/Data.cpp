#include <iostream>
#include <iomanip>
#include <istream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include "Clades.h"
#include "Data.h"
#include "Msg.h"



Data::Data(std::string fileName, std::string errorName, bool useErrors) {

    std::cout << "   * Reading gene presence/absence information" << std::endl;
	/* open the state data file */
	std::ifstream seqStream(fileName.c_str());
	if (!seqStream) 
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	std::string theSequence = "";
	int taxonNum = 0;
	matrix = NULL;
	numTaxa = numChar = 0;
	bool excludeLine = false, charSetLine = false, cladeLine = false;
	bool *tempVec = NULL;
	while( getline(seqStream, linestring).good() )
		{
		std::istringstream linestream(linestring);
		int ch;
		std::string word = "";
		int wordNum = 0;
		int siteNum = 0;
		excludeLine = false;
		charSetLine = false;
        cladeLine   = false;
        Clades* curClade = NULL;
		std::string cmdString = "";
		do
			{
			word = "";
			linestream >> word;
			wordNum++;
			//cout << "word(" << wordNum << ") = " << word << endl;
			if (line == 0)
				{
				/* read the number of taxa/chars from the first line */
				int x;
				std::istringstream buf(word);
				buf >> x;
				if (wordNum == 1)
					numTaxa = x;
				else
					numChar = numSitePatterns = x;
				if (numTaxa > 0 && numChar > 0 && matrix == NULL)
					{	
					matrix = new int*[numTaxa];
					matrix[0] = new int[numTaxa * numChar];
					for (int i=1; i<numTaxa; i++)
						matrix[i] = matrix[i-1] + numChar;
					for (int i=0; i<numTaxa; i++)
						for (int j=0; j<numChar; j++)
							matrix[i][j] = 0;
					isExcluded = new bool[numChar];
					partitionId = new int[numChar];
					patternCount = new int[numChar];
					for (int i=0; i<numChar; i++)
						{
						isExcluded[i] = false;
						partitionId[i] = -1;
						patternCount[i] = 1;
						}
					compressedData = false;
					tempVec = new bool[numChar];
					}
				}
			else
				{
				if (wordNum == 1)
					{
					if ( word == "exclude" )
						excludeLine = true;
					else if ( word == "charset" )
						charSetLine = true;
					else if ( word == "clade" )
						cladeLine = true;
					else
						{
						taxonNames.push_back(word);
						taxonNum++;
						}
					}
				else
					{
                    if (cladeLine == true)
                        {
                        if (wordNum == 2)
                            {
                            curClade = new Clades;
                            clades.push_back(curClade);
                            curClade->setName(word);
                            }
                        else
                            curClade->addTaxon(word);
                        }
                    else
                        {
                        for (int i=0; i<word.length(); i++)
                            {
                            if (siteNum+1 > numChar)
                                Msg::error("Too many characters");
                            char site = word.at(i);
                            matrix[taxonNum-1][siteNum++] = stateID(site);
                            }
                        }
					}
				}
			} while ( (ch=linestream.get()) != EOF );
			
		// NOTE: We probably do not need this bit of code.
		if (line == 0)
			{
			/* the first line should contain the number of taxa and the sequence length */
			std::istringstream buf(word);
			//buf >> genomeSize;
			}
		else
			{
			for (int i=0; i<word.length(); i++)
				{
				char site = word.at(i);
				if (tolower(site) == 'a' || tolower(site) == 'c' || tolower(site) == 'g' || tolower(site) == 't')
					theSequence += tolower(site);
				}
			}
		//cout << linestring << endl;
		line++;
		}	

    if (tempVec != NULL)
        delete [] tempVec;
	
	/* close the file */
	seqStream.close();

    std::cout << "   * Number of user-defined clades: " << clades.size() << std::endl;
    for (std::vector<Clades*>::iterator it = clades.begin(); it != clades.end(); it++)
        (*it)->print();

    // allocate a matrix for the probabilities
    errorProbabilities = new double*[numTaxa];
    errorProbabilities[0] = new double[numTaxa*numChar];
    for (int i=1; i<numTaxa; i++)
        errorProbabilities[i] = errorProbabilities[i-1] + numChar;
    for (int i=0; i<numTaxa; i++)
        for (int j=0; j<numChar; j++)
            errorProbabilities[i][j] = 0.0;

    // read the error file
    if ( useErrors == true )
        {
        std::cout << "   * Reading error probabilities" << std::endl;
        // open the error file
        std::ifstream errStream(errorName.c_str());
        if (!errStream)
            {
            std::cerr << "Cannot open file \"" + errorName + "\"" << std::endl;
            exit(1);
            }

        linestring = "";
        line = 0;
        taxonNum = 0;
        while( getline(errStream, linestring).good() )
            {
            std::istringstream linestream(linestring);
            int ch;
            std::string word = "";
            int wordNum = 0;
            int siteNum = 0;
            excludeLine = false;
            charSetLine = false;
            std::string cmdString = "";
            do
                {
                word = "";
                linestream >> word;
                wordNum++;
                //std::cout << "word(" << wordNum << ") = " << word << std::endl;
                std::string numStr = "";
                if (line == 0)
                    {
                    /* read the number of taxa/chars from the first line */
                    int x;
                    std::istringstream buf(word);
                    buf >> x;
                    if (wordNum == 1)
                        {
                        if (x != numTaxa)
                            Msg::error("Incorrect number of taxa in error file");
                        }
                    else
                        {
                        if (x != numChar)
                            Msg::error("Incorrect number of characters in error file");
                        }
                    }
                else
                    {
                    if (wordNum == 1)
                        {
                        if (word != taxonNames[taxonNum])
                            Msg::error("Mismatched taxon names in error file");
                        taxonNum++;
                        std::cout << "     Reading taxon " << taxonNum << std::endl;
                        }
                    else
                        {
                        for (int i=0; i<word.length(); i++)
                            {
                            if (word.at(i) != ' ' && word.at(i) != ',')
                                {
                                numStr += word.at(i);
                                }
                            else
                                {
                                if (numStr != "")
                                    {
                                    double x;
                                    std::istringstream buf(word);
                                    buf >> x;
                                    errorProbabilities[taxonNum-1][siteNum++] = x;
                                    numStr = "";
                                    }
                                }
                            }
                        }
                    }
                    
                // check that the numStr is empty as it should be at this point (I think)
                if (numStr != "")
                    {
                    double x;
                    std::istringstream buf(word);
                    buf >> x;
                    errorProbabilities[taxonNum-1][siteNum++] = x;
                    numStr = "";
                    }

                } while ( (ch=linestream.get()) != EOF );
                
            line++;
            }	
        
        // close the error file
        errStream.close();
        }
    
    // check the data matrix for all zero or all one patterns, which should not be there
    std::cout <<"   * Checking data for illegal patterns" << std::endl;
    bool* validPattern = new bool[numChar];
    int numInvalidPatterns = 0;
    for (int c=0; c<numChar; c++)
        {
        validPattern[c] = true;
        int numZeros = 0, numOnes = 0;
        for (int i=0; i<numTaxa; i++)
            {
            int charCode = getCharacter(i, c);
            int possibleChars[2] = { 0, 0 };
            getPossibleChars(charCode, possibleChars);
            numZeros += possibleChars[0];
            numOnes  += possibleChars[1];
            }
#       if defined(CONDITION_NO_CONSTANT)
        if ( (numZeros == 0 || numZeros == numTaxa) || (numOnes == 0 || numOnes == numTaxa) )
            {
            /*std::string errStr = "Removing illegal (all 0 or all 1) character patterns in data matrix";
            std::cout << c << " -- ";
            std::vector<int> x = getSitePattern(c);
            for (int z=0; z<x.size(); z++)
                std::cout << x[z];
            std::cout << std::endl;*/
            validPattern[c] = false;
            numInvalidPatterns++;
            }
#       elif defined(CONDITION_NO_ZEROS)
        if ( numZeros == numTaxa )
            {
            /*std::string errStr = "Removing illegal (all 0) character patterns in data matrix";
            std::cout << c << " -- ";
            std::vector<int> x = getSitePattern(c);
            for (int z=0; z<x.size(); z++)
                std::cout << x[z];
            std::cout << std::endl;*/
            validPattern[c] = false;
            numInvalidPatterns++;
            }
#       endif
            
        }
    
    if (numInvalidPatterns == 0)
        std::cout <<"     No invalid patterns found" << std::endl;
    else if (numInvalidPatterns == 1)
        std::cout <<"     Discarded one invalid pattern" << std::endl;
    else
        std::cout <<"     Discarded " << numInvalidPatterns << " invalid patterns" << std::endl;
    
    int** newIntMatrix       = allocateIntDataMatrix(numTaxa, numChar-numInvalidPatterns);
    double** newDoubleMatrix = allocateDoubleDataMatrix(numTaxa, numChar-numInvalidPatterns);
    for (int j=0, k=0; j<numChar; j++)
        {
        if (validPattern[j] == true)
            {
            for (int i=0; i<numTaxa; i++)
                {
                newIntMatrix[i][k] = matrix[i][j];
                newDoubleMatrix[i][k] = errorProbabilities[i][j];
                }
            k++;
            }
        }
    freeDataMatrix(matrix);
    freeDataMatrix(errorProbabilities);
    matrix = newIntMatrix;
    errorProbabilities = newDoubleMatrix;
    numChar -= numInvalidPatterns;

}



Data::~Data(void) {

	delete [] matrix[0];
	delete [] matrix;
	delete [] patternCount;
	delete [] isExcluded;
	delete [] partitionId;
	if (compressedData == true)
		{
		delete [] compressedMatrix[0];
		delete [] compressedMatrix;
		}
	delete [] errorProbabilities[0];
	delete [] errorProbabilities;
}


int** Data::allocateIntDataMatrix(int nt, int nc) {

    int** m;
    m = new int*[nt];
    m[0] = new int[nt * nc];
    for (int i=1; i<nt; i++)
        m[i] = m[i-1] + nc;
    for (int i=0; i<nt; i++)
        for (int j=0; j<nc; j++)
            m[i][j] = 0;
    return m;
}

double** Data::allocateDoubleDataMatrix(int nt, int nc) {

    double** m;
    m = new double*[nt];
    m[0] = new double[nt * nc];
    for (int i=1; i<nt; i++)
        m[i] = m[i-1] + nc;
    for (int i=0; i<nt; i++)
        for (int j=0; j<nc; j++)
            m[i][j] = 0;
    return m;
}

void Data::freeDataMatrix(int** m) {

    delete [] m[0];
    delete [] m;
}

void Data::freeDataMatrix(double** m) {

    delete [] m[0];
    delete [] m;
}

void Data::compress(void) {

	if (compressedData == false)
		{
		int *tempCount = new int[numChar];
		for (int j=0; j<numChar; j++)
			tempCount[j] = 1;
		
		for (int j=0; j<numChar; j++)
			{
			if (tempCount[j] > 0)
				{
				for (int k=j+1; k<numChar; k++)
					{
					if (tempCount[k] == 0)
						continue;
					bool isSame = true;
					for (int i=0; i<numTaxa; i++)
						{
						if (matrix[i][j] != matrix[i][k])
							{
							isSame = false;
							break;
							}
						}
					if (isSame == true)
						{
						tempCount[j]++;
						tempCount[k] = 0;
						}
					}
				}
			}
		numSitePatterns = 0;
		for (int j=0; j<numChar; j++)
			if (tempCount[j] > 0)
				numSitePatterns++;
				
		compressedMatrix = new int*[numTaxa];
		compressedMatrix[0] = new int[numTaxa * numSitePatterns];
		for (int i=1; i<numTaxa; i++)
			compressedMatrix[i] = compressedMatrix[i-1] + numSitePatterns;
		for (int i=0; i<numTaxa; i++)
			for (int j=0; j<numSitePatterns; j++)
				compressedMatrix[i][j] = 0;
		
		delete [] patternCount;
		patternCount = new int[numSitePatterns];
		for (int j=0, k=0; j<numChar; j++)
			{
			if (tempCount[j] > 0)
				{
				for (int i=0; i<numTaxa; i++)
					compressedMatrix[i][k] = matrix[i][j];
				patternCount[k] = tempCount[j];
				k++;
				}
			}
		compressedData = true;
		delete [] tempCount;
		}
		
}



int Data::getCharacter(int i, int j) {

	if (compressedData == true)
		return compressedMatrix[i][j];
	else
		return matrix[i][j];
		
}

double Data::getErrorProbability (int i, int j) {

    return errorProbabilities[i][j];
}



/*-------------------------------------------------------------------
|
|   GetPossibleNucs: 
|
|   This function initializes a vector, nuc[MAX_NUM_STATES]. The four elements
|   of nuc correspond to the four nucleotides in alphabetical order.
|   We are assuming that the nucCode is a binary representation of
|   the nucleotides that are consistent with the observation. For
|   example, if we observe an A, then the nucCode is 1 and the 
|   function initalizes nuc[0] = 1 and the other elements of nuc
|   to be 0.
|
|   Observation    nucCode        nuc
|        A            1           1000
|        C            2           0100
|        G            4           0010
|        T            8           0001
|        R            5           1010
|        Y           10           0101
|        M            3           1100
|        K           12           0011
|        S            6           0110
|        W            9           1001
|        H           11           1101
|        B           14           0111
|        V            7           1110
|        D           13           1011
|        N - ?       15           1111
|
-------------------------------------------------------------------*/
void Data::getPossibleChars(int charCode, int* possibleChars) {

	if (charCode == 1)
		{
		possibleChars[0] = 1;
		possibleChars[1] = 0;
		}
	else if (charCode == 2)
		{
		possibleChars[0] = 0;
		possibleChars[1] = 1;
		}
	else if (charCode == 3)
		{
		possibleChars[0] = 1;
		possibleChars[1] = 1;
		}
}

void Data::listTaxa(void) {

	int i = 1;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		std::cout << std::setw(4) << i++ << " -- " << (*p) << std::endl;

}

int Data::getNumSubsets(void) {

	int largestId = 0;
	for (int i=0; i<numChar; i++)
		{
		if (partitionId[i] > largestId)
			{
			largestId = partitionId[i];
			}
		}
		
	bool *isPartHere = new bool[largestId];
	for (int i=0; i<largestId; i++)
		isPartHere[i] = false;
	for (int i=0; i<numChar; i++)
		{
		if (isExcluded[i] == false)
			{
			isPartHere[ partitionId[i]-1 ] = true;
			}
		}
	int numParts = 0;
	for (int i=0; i<largestId; i++)
		if (isPartHere[i] == true)
			numParts++;
	delete [] isPartHere;
		
	return numParts;
		
}

std::vector<int> Data::getSitePattern(int idx) {

    std::vector<int> x;
    for (int i=0; i<numTaxa; i++)
        {
        x.push_back( matrix[i][idx] );
        }
    return x;
}

std::string Data::getTaxonName(int i) {

	return taxonNames[i];

}

int Data::getTaxonIndex(std::string ns) {

	int taxonIndex = -1;
	int i = 0;
	for (std::vector<std::string>::iterator p=taxonNames.begin(); p != taxonNames.end(); p++)
		{
		if ( (*p) == ns )
			{
			taxonIndex = i;
			break;
			}
		i++;
		}
	return taxonIndex;

}

void Data::interpretString(std::string s, bool *v, int n) {

	for (int i=0; i<n; i++)
		v[i] = false;
	int nums[3];
	int numToRemember = 0;
	//cout << "s = \"" << s << "\"" << endl;

	/* push the individual words (numbers, hyphens, or back slashes into a vector */
	std::vector<std::string> cmds;
	std::istringstream linestream(s);
	int ch;
	std::string word = "";
	do
		{
		word = "";
		linestream >> word;
		if (word != "")
			{
			cmds.push_back( word );
			}
		} while( (ch=linestream.get()) != EOF );
		
	/* loop over the vector of individual words and correctly interpret things */
	int i = 0;
	for (std::vector<std::string>::iterator p=cmds.begin(); p != cmds.end(); p++)
		{
		//cout << "\"" << (*p) << "\"" << endl;
		if ( isdigit((*p)[0]) )
			{
			
			/* we can potentially complete a sentence */
			std::string prevWord = "";
			if (i > 0)
				prevWord = cmds[i-1];
			std::string nextWord = "";
			if (i != cmds.size() - 1)
				nextWord = cmds[i+1];
			//cout << "\"" << prevWord << "\" ** \"" << cmds[i] << "\" ** \"" << nextWord << "\"" << endl;
			
			int x;
			std::istringstream buf(cmds[i]);
			buf >> x;

			if (prevWord == "" || isNumber(prevWord) == true)
				{
				nums[0] = x;
				numToRemember = 1;
				}
			else if (prevWord == "-")
				{
				nums[1] = x;
				numToRemember = 2;
				}
			else if (prevWord == "\\")
				{
				nums[2] = x;
				numToRemember = 3;
				}
			else
				{
				std::cerr << "ERROR: Problem interpreting string" << std::endl;
				exit(1);
				}
			
			if ( (prevWord == "" || isNumber(prevWord) == true) && (nextWord == "" || isNumber(nextWord) == true) )
				{
				v[nums[0]-1] = true;
				numToRemember = 0;
				}
			else if ( prevWord == "-" && (nextWord == "" || isNumber(nextWord) == true) )
				{
				for (int i=nums[0]-1; i<nums[1]; i++)
					v[i] = true;
				numToRemember = 0;
				}
			else if ( prevWord == "\\" )
				{
				for (int i=nums[0]-1, k=nums[2]; i<nums[1]; i++, k++)
					if ( k % nums[2] == 0 )
						v[i] = true;
				numToRemember = 0;
				}
			}
		i++;
		}
#	if 0
	for (int i=0; i<n; i++)
		{
		if (v[i] == true)
			cout << i+1 << " ";
		}
	cout << endl;
#	endif

}

bool Data::isNumber(std::string s) {

	if (s == "")
		return false;

	bool isnum = true;
	for (int i=0; i<s.size(); i++)
		if ( !isdigit(s[i]) )
			isnum = false;
	return isnum;
	
}

int Data::stateID(char st) {
	
	if (st == '0')
		{
		return 1;
		}
	else if (st == '1')
		{
		return 2;
		}
	else if (st == '?' || st == 'N')
		{
		return 3;
		}
	else
		return -1;
}

void Data::print(void) {

	int **x;
	if (compressedData == false)
		x = matrix;
	else
		x = compressedMatrix;
		
	for (int j=0; j<numSitePatterns; j++)
		{
		std::cout << std::setw(4) << j << " -- ";
		for (int i=0; i<numTaxa; i++)
			{
			std::cout << x[i][j];
			}
		std::cout << std::endl;
		}
		
	int sum = 0;
	for (int j=0; j<numSitePatterns; j++)
		sum += patternCount[j];
	std::cout << "Number of sites = " << sum << std::endl;
}

void Data::printErrorProbabilities(void) {

    for (int j=0; j<numChar; j++)
        {
        std::cout << std::setw(4) << j << " -- ";
        for (int i=0; i<numTaxa; i++)
            {
            std::cout << std::fixed << std::setprecision(3) << errorProbabilities[i][j] << " ";
            }
        std::cout << std::endl;
        }
}

void Data::uncompress(void) {

	if (compressedData == true)
		{
		delete [] compressedMatrix[0];
		delete [] compressedMatrix;
		delete [] patternCount;
		patternCount = new int[numChar];
		for (int j=0; j<numChar; j++)
			patternCount[j] = 1;
		numSitePatterns = numChar;
		compressedData = false;
		}

}

