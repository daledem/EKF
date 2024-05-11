#include "../include/TimeUpdate.h"

Matrix TimeUpdate(const Matrix &P,const Matrix &Phi, double Qdt) {
    return Phi*P*Phi.trans() + Qdt;
}
