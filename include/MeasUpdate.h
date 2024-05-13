#ifndef MEASUPDATE_H
#define MEASUPDATE_H

#include "./Matrix.h"

void MeasUpdate(Matrix& K,Matrix& xOut,Matrix& POut,Matrix& x, double z,double g,double s,Matrix& G, Matrix& P,int n);

#endif //MEASUPDATE_H
