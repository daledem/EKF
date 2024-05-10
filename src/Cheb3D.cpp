#include "Cheb3D.h"

Matrix Cheb3D(double t, int N, double Ta, double Tb,const Matrix &Cx,const Matrix &Cy,const Matrix &Cz) {
    // Check validity
    if ( (t<Ta) || (Tb<t) )
        throw "ERROR: Time out of range in Cheb3D::Value\n";

    // Clenshaw algorithm
    double tau = (2*t-Ta-Tb)/(Tb-Ta);

    Matrix f1(1,3);
    Matrix f2(1,3);
    Matrix old_f1(1,3);

    Matrix Caux(1,3);

    for(int i = N; i >= 2; --i){
        old_f1 = f1;

        Caux(1,1) = Cx(1,i);
        Caux(1,2) = Cy(1,i);
        Caux(1,3) = Cz(1,i);

        f1 = f1*2*tau-f2+Caux;
        f2 = old_f1;
    }

    Caux(1,1) = Cx(1,1);
    Caux(1,2) = Cy(1,1);
    Caux(1,3) = Cz(1,1);

    return f1*tau-f2+Caux;
}
