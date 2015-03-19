
#define CONDITION_NO_ZEROS
#undef CONDITION_NO_CONSTANT

#ifndef Data_H
#define Data_H

#define PRESENCE_PROBABILITY_DELIMITER ','

#include <string>
#include <vector>

class Clades;

class Data {

	public:
                                    Data(std::string fileName, std::string genePresenceProbabilitiesFileName, bool useGenePresenceProbabilities);
                                   ~Data(void);
        void                        compress(void);
        int                         getNumTaxa(void) { return numTaxa; }
        int                         getNumChar(void) { return (compressedData == true ? numSitePatterns : numChar); }
        void                        getPossibleChars(int charCode, int* possibleChars);
        int                         getCharacter(int i, int j);
        double                      getGenePresenceProbability (int i, int j);
        int                         getTaxonIndex(std::string ns);
        int                         getPartitionId(int i) { return partitionId[i]; }
        void                        listTaxa(void);
        void                        print(void);
        void                        printGenePresenceProbabilities(void);
        void                        uncompress(void);
        Clades*                     getClade(int idx) { return clades[idx]; }
        std::string                 getTaxonName(int i);
        int                         getNumOfPattern(int i) { return patternCount[i]; }
        bool                        getIsExcluded(int i) { return isExcluded[i]; }
        int                         getNumClades(void) { return (int)clades.size(); }
        int                         getNumSubsets(void);
        std::vector<int>            getSitePattern(int idx);

    private:
        int**                       allocateIntDataMatrix(int nt, int nc);
        double**                    allocateDoubleDataMatrix(int nt, int nc);
        void                        freeDataMatrix(int** m);
        void                        freeDataMatrix(double** m);
        void                        interpretString(std::string s, bool *v, int n);
        bool                        isNumber(std::string s);
        int                         stateID(char st);
        int                         numTaxa;
        int                         numChar;
        int                         numSitePatterns;
        std::vector<std::string>    taxonNames;
        bool                        compressedData;
        bool*                       isExcluded;
        int*                        partitionId;
        int**                       matrix;
        int**                       compressedMatrix;
        int*                        patternCount;
        double**                    genePresenceProbabilities;
        std::vector<Clades*>        clades;
};

#endif
