//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/13
//
//------------------------------------------------------------------------------
#include "../include/MeasUpdate.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void MeasUpdate(Matrix &K, Matrix& x, double z, double g, double s,const Matrix &G, Matrix &P, int n)
//------------------------------------------------------------------------------
/**
 *   I don´t know what this does
 *
 * @param[out] <K> Kalman gain
 * @param[in,out] <x> State update
 * @param <z> double
 * @param <g> double
 * @param <s> double
 * @param <G> Matrix
 * @param[in,out] <P> Covariance update
 * @param <n> int
 *
 */
//------------------------------------------------------------------------------

void MeasUpdate(Matrix &K, Matrix& x, double z, double g, double s,const Matrix &G, Matrix &P, int n) {
    int m;
    m = 1;
    Matrix Inv_W(m,m);

    for(int i = 1; i<= m;i++)
        Inv_W(i,i) = s*s;    // Inverse weight (measurement covariance)

    // Kalman gain
    K = P*G.trans()*(Inv_W+G*P*G.trans()).inverse();

    // State update
    x = x + K*(z-g);

    // Covariance update
    P = (Matrix::eye(n) - ((K*G)))*P;
}

