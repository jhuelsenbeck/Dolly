#include <iostream>
#include "Clades.h"



Clades::Clades(void) {

}

Clades::~Clades(void) {

}

void Clades::print(void) {

    std::cout << "     Clade: " << name << std::endl;
    for (int i=0; i<taxaInClade.size(); i++)
        std::cout << "        " << taxaInClade[i] << std::endl;
}