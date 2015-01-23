#ifndef Clades_H
#define Clades_H

#include <string>
#include <vector>

class Clades {

	public:
                                    Clades(void);
                                   ~Clades(void);
        void                        addTaxon(std::string s) { taxaInClade.push_back(s); }
        std::string                 getName(void) { return name; }
        std::string                 getTaxon(int idx) { return taxaInClade[idx]; }
        void                        print(void);
        void                        setName(std::string s) { name = s; }

    private:
        std::string                 name;
        std::vector<std::string>    taxaInClade;
};


#endif
