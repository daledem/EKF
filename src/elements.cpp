#include "../include/elements.h"

void elements(double& p, double& a, double& e, double& i, double& Omega, double& omega, double& M, const Matrix& y){
    double pi2,magh,H,u,R,eCosE,eSinE,e2,E,nu;

    pi2 = 2*Const::pi;

    Matrix r = y.getColumnaByIndex(1).trans();                                        // Position
    Matrix v = y.getColumnaByIndex(2).trans();                                        // Velocity

    Matrix h = Matrix::cross(r,v);                                    // Areal velocity
    magh = Matrix::norm(h);
    p = magh*magh/Const::GM_Earth;
    H = Matrix::norm(h);

    Omega = atan2 ( h(1,1), -h(2) );                     // Long. ascend. node
    Omega = fmod(Omega,pi2);
    i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); // Inclination
    u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    // Arg. of latitude

    R  = Matrix::norm(r);                                      // Distance

    a = 1/(2/R-Matrix::dot(v,v)/Const::GM_Earth);               // Semi-major axis

            eCosE = 1-R/a;                                     // e*cos(E)
    eSinE = Matrix::dot(r,v)/sqrt(Const::GM_Earth*a);           // e*sin(E)

    e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                     // Eccentricity
    E  = atan2(eSinE,eCosE);                           // Eccentric anomaly

    M  = fmod(E-eSinE,pi2);                             // Mean anomaly

    nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);          // True anomaly

    omega = fmod(u-nu,pi2);                             // Arg. of perihelion
}