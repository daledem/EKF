#ifndef DEINTEG_H
#define DEINTEG_H

#include <cmath>
#include <algorithm>
#include "./Matrix.h"
#include "./sign_.h"
#include "./SAT_Const.h"

Matrix DEInteg(Matrix (*func)(double x,const Matrix& Y),double t,double tout,double relerr,double abserr,int n_eqn,Matrix y);

#endif //DEINTEG_H
