#ifndef CharacterChanges_H
#define CharacterChanges_H

#define DIMENSIONS 40

#include <string>

class CharacterChanges {

	public:
                                    CharacterChanges(int mc);
        std::string                 getChangeProbsHeader(void);
        std::string                 getChangeProbs(void);
        void                        updateChanges(int n01, int n10);
        double                      averageNumberGains(void);
        double                      averageNumberLosses(void);
        void                        print(void);

    private:
        int                         dim;
        int                         n;
        int                         myChar;
        int                         changes[DIMENSIONS][DIMENSIONS];
};


#endif
