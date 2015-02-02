#include <iomanip>
#include <iostream>
#include <stdio.h>
#include "CharacterChanges.h"
#include "Msg.h"



CharacterChanges::CharacterChanges(int mc) {

    myChar = mc;
    dim = DIMENSIONS;
    n = 0;
    for (int i=0; i<dim; i++)
        for (int j=0; j<dim; j++)
            changes[i][j] = 0.0;
}

double CharacterChanges::averageNumberGains(void) {

    int sum = 0;
    for (int i=0; i<dim; i++)
        {
        for (int j=0; j<dim; j++)
            sum += i * changes[i][j];
        }
    return (double)sum / n;
}

double CharacterChanges::averageNumberLosses(void) {

    int sum = 0;
    for (int j=0; j<dim; j++)
        {
        for (int i=0; i<dim; i++)
            sum += j * changes[i][j];
        }
    return (double)sum / n;
}

std::string CharacterChanges::getChangeProbsHeader(void) {

    std::string str = "";
    for (int i=0; i<dim; i++)
        {
        for (int j=0; j<dim; j++)
            {
            char cStr[20];
            sprintf(cStr, "P(%d,%d)\t", i, j);
            str += cStr;
            }
        }
    return str;
}

std::string CharacterChanges::getChangeProbs(void) {

    std::string str = "";
    for (int i=0; i<dim; i++)
        {
        for (int j=0; j<dim; j++)
            {
            char cStr[20];
            sprintf(cStr, "%1.3lf\t", (double)changes[i][j] / n);
            str += cStr;
            }
        }
    return str;
}

void CharacterChanges::updateChanges(int n01, int n10) {

    if (n01 > dim-1 || n10 > dim-1)
        {
        Msg::warning("Too many changes to record");
        std::cout << myChar << " -- " << n01 << " " << n10 << std::endl;
        return;
        }
    changes[n01][n10]++;
    n++;
}
