//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/11
//
//------------------------------------------------------------------------------
#include "../include/VarEqn.h"

extern auxParam AuxParam;
extern Matrix eopdata;

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix VarEqn(double x, const Matrix &yPhi)
//------------------------------------------------------------------------------
/**
 *   Computes the variational equations, i.e. the derivative of the state vector
 *      and the state transition matrix
 *
 * @param <x> Time since epoch in [s]
 * @param <yPhi> (6+36)-dim vector comprising the state vector (y) and the
 *                  state transition matrix (Phi) in column wise storage order
 *
 * @return Derivative of yPhi
 *
 */
//------------------------------------------------------------------------------
Matrix VarEqn(double x, const Matrix &yPhi) {
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,
    UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,Mjd_UT1;
    IERS(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,eopdata,AuxParam.Mjd_UTC,'l');
    timediff(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC,UT1_UTC,TAI_UTC);
    Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    // Transformation matrix
    Matrix P = PrecMatrix(Const::MJD_J2000,AuxParam.Mjd_TT + x/86400);
    Matrix N = NutMatrix(AuxParam.Mjd_TT + x/86400);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r = yPhi.getColumnaByIndex(1,1,3);
    Matrix v = yPhi.getColumnaByIndex(1,4,6);
    Matrix Phi(6,6);

    // State transition matrix
    for (int j = 1; j <= 6; j++)
        for(int i = 1; i <= 6 ; i++)
            Phi(i,j) = yPhi(6*j+i,1);


    // Acceleration and gradient
    Matrix a = AccelHarmonic ( r, E, AuxParam.n, AuxParam.m );
    Matrix G = G_AccelHarmonic ( r, E, AuxParam.n, AuxParam.m );

    // Time derivative of state transition matrix
    Matrix yPhip(42,1);
    Matrix dfdy(6,6);

    for (int i = 1; i<= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i,j) = 0.0;                 // dv/dr(i,j)
            dfdy(i+3,j) = G(i,j);            // da/dr(i,j)
            if ( i==j )
                dfdy(i,j+3) = 1;
            else
                dfdy(i,j+3) = 0;             // dv/dv(i,j)

            dfdy(i+3,j+3) = 0.0;             // da/dv(i,j)
        }
    }

    Matrix Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix
    for (int i = 1; i <= 3 ; i++) {
        yPhip(i,1)   = v(i,1);                 // dr/dt(i)
        yPhip(i+3,1) = a(i,1);                 // dv/dt(i)
    }

    for (int i = 1; i <= 6; i++)
        for (int j = 1; j <= 6; j++)
            yPhip(6*j+i,1) = Phip(i,j);     // dPhi/dt(i,j)

    return yPhip;
}
