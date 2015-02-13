#ifndef Settings_H
#define Settings_H

#include <string>


class Settings {

	public:
                        Settings(int argc, char* argv[]);
        int             getBurnIn(void) { return burnIn; }
        int             getChainLength(void) { return chainLength; }
        std::string     getGenePresenceProbabilitiesFileName(void) { return genePresenceProbabilitiesFileName; }
        bool            getUseGenePresenceProbs(void) { return useGenePresenceProbs; }
        std::string     getInputFileName(void) { return inputFileName; }
        std::string     getOutputFileName(void) { return outputFileName; }
        int             getPrintFrequency(void) { return printFrequency; }
        int             getSampleFrequency(void) { return sampleFrequency; }
        std::string     getTreeFileName(void) { return treeFileName; }
        double          getBranchLengthPrior(void) { return branchLengthPrior; }
        bool            getFixBranchProportionsToUserTree(void) { return fixBranchProportionsToUserTree; }
        void            print(void);
        void            setBurnIn(int x) { burnIn = x; }
        void            setChainLength(int x) { chainLength = x; }
        void            setInputFileName(std::string s) { inputFileName = s; }
        void            setOutputFileName(std::string s) { outputFileName = s; }
        void            setPrintFrequency(int x) { printFrequency = x; }
        void            setSampleFrequency(int x) { sampleFrequency = x; }
        void            setTreeFileName(std::string s) { treeFileName = s; }
        void            setGenePresenceProbabilitiesFileName(std::string s) { genePresenceProbabilitiesFileName = s; }
        void            setUseGenePresenceProbs(bool tf) { useGenePresenceProbs = tf; }
        void            setBranchLengthPrior(double x) { branchLengthPrior = x; }
        void            setBranchProportionsToUserTree(bool tf) { fixBranchProportionsToUserTree = tf; }
    private:
        void            printUsage(void);
        int             chainLength;
        int             burnIn;
        bool            useGenePresenceProbs;
        std::string     inputFileName;
        std::string     genePresenceProbabilitiesFileName;
        std::string     treeFileName;
        std::string     outputFileName;
        int             printFrequency;
        int             sampleFrequency;
        double          branchLengthPrior;
        bool            fixBranchProportionsToUserTree;
};


#endif
