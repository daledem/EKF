#ifndef MEASUPDATE_H
#define MEASUPDATE_H

#include "./Matrix.h"

void MeasUpdate(Matrix& K,Matrix& xOut,Matrix& POut,const Matrix& x, double z,double g,double s,const Matrix& G,const Matrix& P,int n);

#endif //MEASUPDATE_H
