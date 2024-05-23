#include "../include/angl.h"

double angl(const Matrix &vec1, const Matrix &vec2) {
    double small, undefined, magv1,magv2,temp;

    small     = 0.00000001;
    undefined = 999999.1;

    magv1 = Matrix::norm(vec1);
    magv2 = Matrix::norm(vec2);

    if (magv1*magv2 > small^2) {
        temp= Matrix::dot(vec1,vec2) / (magv1*magv2);
        if (fabs( temp ) > 1.0) {
            temp= sign(temp) * 1.0;
        }
        return acos( temp );
    }

    return undefined;
}
