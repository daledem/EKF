//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/23
//
//------------------------------------------------------------------------------
#include "../include/elements.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// void elements(double& p, double& a, double& e, double& i, double& Omega,
//                  double& omega, double& M, const Matrix& y)
//------------------------------------------------------------------------------
/**
 *   Computes the osculating Keplerian elements from the satellite state
 *   vector for elliptic orbits
 *
 * @param[out] <p> semilatus rectum [m]
 * @param[out] <a> Semimajor axis
 * @param[out] <e> Eccentricity
 * @param[out] <i> Inclination [rad]
 * @param[out] <Omega> Longitude of the ascending node [rad]
 * @param[out] <omega> Argument of pericenter [rad]
 * @param[out] <M> Mean anomaly [rad]
 * @param <y> State vector (x,y,z,vx,vy,vz)
 *
 * @note The function cannot be used with state vectors describing a circular
 *         or non-inclined orbit.
 */
//------------------------------------------------------------------------------
void elements(double& p, double& a, double& e, double& i, double& Omega, double& omega, double& M, const Matrix& y){
    double pi2,magh,H,u,R,eCosE,eSinE,e2,E,nu;

    pi2 = 2*Const::pi;

    Matrix r = y.getColumnaByIndex(1).trans();                                        // Position
    Matrix v = y.getColumnaByIndex(2).trans();                                        // Velocity

    Matrix h = Matrix::cross(r,v);                                    // Areal velocity
    magh = Matrix::norm(h);
    p = magh*magh/Const::GM_Earth;
    H = Matrix::norm(h);

    Omega = atan2 ( h(1,1), -h(1,2) );                     // Long. ascend. node
    Omega = mod(Omega,pi2);
    i     = atan2 ( sqrt(h(1,1)*h(1,1)+h(1,2)*h(1,2)), h(1,3) ); // Inclination
    u     = atan2 ( r(1,3)*H, -r(1,1)*h(1,2)+r(1,2)*h(1,1) );    // Arg. of latitude

    R  = Matrix::norm(r);                                      // Distance

    a = 1/(2/R-Matrix::dot(v,v)/Const::GM_Earth);               // Semi-major axis

    eCosE = 1-R/a;                                     // e*cos(E)
    eSinE = Matrix::dot(r,v)/sqrt(Const::GM_Earth*a);           // e*sin(E)

    e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                     // Eccentricity
    E  = atan2(eSinE,eCosE);                           // Eccentric anomaly

    M  = mod(E-eSinE,pi2);                             // Mean anomaly

    nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          // True anomaly

    omega = mod(u-nu,pi2);                             // Arg. of perihelion
}