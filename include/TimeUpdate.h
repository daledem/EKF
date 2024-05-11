#ifndef TIMEUPDATE_H
#define TIMEUPDATE_H

#include "./Matrix.h"

Matrix TimeUpdate(const Matrix& P,const Matrix& Phi, double Qdt = 0.0);

#endif //TIMEUPDATE_H
