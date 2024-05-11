#include "../include/sign_.h"
#include <cmath>

double sign_(double a,double b){
    double result;
    if (b>=0.0)
        result = fabs(a);
    else
        result = - fabs(a);
    return result;
}
