#include "../include/mod.h"

// Codigo obtenido de https://stackoverflow.com/questions/28888619/modulo-function-in-c-that-behaves-like-mod-in-matlab
double mod(double a, double b) {
    double result = fmod(a, b);
    return result >= 0 ? result : result + b;
}
