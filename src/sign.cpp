#include "../include/sign.h"

double sign(double a) {
    if(a != 0)
        return a/fabs(a);
    return 0;
}

