#ifndef PROYECTO_UNIT_H
#define PROYECTO_UNIT_H

#include "./Matrix.h"

//------------------------------------------------------------------------------
// unit(const Matrix& vec)
//------------------------------------------------------------------------------
/**
*   this function calculates a unit vector given the original vector. if a
*       zero vector is input, the vector is set to zero.
*
* @param <vec> vector
*
* @return unit vector
*
*/
//------------------------------------------------------------------------------

Matrix unit(const Matrix& vec);

#endif //PROYECTO_UNIT_H
