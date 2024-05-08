#ifndef PROYECTO_GLOBALS_H
#define PROYECTO_GLOBALS_H

#include "Matrix.h"
#include <cstdio>
#include <cstdlib>

class Global{


public:
    static Matrix *Cnm;
    static Matrix *Snm;

    static Matrix *eopdata;

    static struct {
        double Mjd_UTC;
        int n;
        int m;
        int sun;
        int moon;
        int planets;
        double Mjd_TT;
    }AuxParam;

    static int n_eqn;

    static void GGM03S();
    static void eop19620101(int c = 21413);
};

#endif //PROYECTO_GLOBALS_H
