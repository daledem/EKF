#include "../include/GHAMatrix.h"

Matrix GHAMatrix(double Mjd_UT1) {
    return R_z( gast(Mjd_UT1) );
}
