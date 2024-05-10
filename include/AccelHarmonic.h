#ifndef PROYECTO_ACCELHARMONIC_H
#define PROYECTO_ACCELHARMONIC_H

#include <cmath>
#include "./Matrix.h"
#include "./Legendre.h"
#include "./SAT_Const.h"

Matrix AccelHarmonic(const Matrix& r,const Matrix& E,int n_max,int m_max);

#endif //PROYECTO_ACCELHARMONIC_H
