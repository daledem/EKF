#ifndef CHEB3D_H
#define CHEB3D_H

#include "Matrix.h"

Matrix Cheb3D(double t,int N,double Ta,double Tb,const Matrix& Cx,const Matrix& Cy,const Matrix& Cz);

#endif //CHEB3D_H
