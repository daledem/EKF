#include "../include/LTC.h"

Matrix LTC(double lon, double lat) {
    double Aux;
    Matrix M = R_y(-1.0*lat)*R_z(lon);

    for(int j = 1; j <= 3; j++) {
        Aux=M(1,j);
        M(1,j)=M(2,j);
        M(2,j)=M(3,j);
        M(3,j)= Aux;
    }

    return M;
}
