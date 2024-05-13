#include "../include/MeasUpdate.h"

void MeasUpdate(Matrix &K, Matrix &xOut, Matrix &POut,Matrix& x, double z, double g, double s, Matrix &G, Matrix &P, int n) {
    int m;
    m = 1;
    Matrix Inv_W(m,m);

    for(int i = 1; i<= m;i++)
        Inv_W(i,i) = s*s;    // Inverse weight (measurement covariance)

    // Kalman gain
    K = P*G.trans()*(Inv_W+G*P*G.trans()).inverse();

    // State update
    xOut = x + K*(z-g);

    // Covariance update
    POut = (Matrix::eye(n)-K*G)*P;
}

