#include "../include/AccelPointMass.h"

Matrix AccelPointMass(Matrix& r,Matrix& s, double GM) {

    Matrix d = r - s;

    // Acceleration 
    return -GM * ( d/(pow(Matrix::norm(d),3)) + s/(pow(Matrix::norm(s),3)) );
}


