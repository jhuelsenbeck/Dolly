#ifndef Mcmc_H
#define Mcmc_H

#include <fstream>

class MbRandom;
class Model;
class Settings;

class Mcmc {

	public:
                                                Mcmc(Settings* sp, Model* mp, MbRandom* rp);
                                               ~Mcmc(void);

    private:
        void                                    runChain(void);
        void                                    saveStates(int n, double lnL);
        Model*                                  modelPtr;
        MbRandom*                               ranPtr;
        Settings*                               settingsPtr;
        int                                     chainLength;
        int                                     burnIn;
        int                                     printFrequency;
        int                                     sampleFrequency;
        std::ofstream                           treeOut;
        std::ofstream                           parmOut;
};


#endif
