#ifndef PROYECTO_GLOBALS_H
#define PROYECTO_GLOBALS_H

#include "Matrix.h"
#include <cstdio>
#include <cstdlib>

class Global{
public:
    static Matrix *eopdata;

    static void eop19620101(int c);
};

#endif //PROYECTO_GLOBALS_H
