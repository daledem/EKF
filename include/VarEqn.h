#ifndef VAREQN_H
#define VAREQN_H

#include "./Matrix.h"
#include "./SAT_Const.h"
#include "./Globals.h"
#include "./IERS.h"
#include "./timediff.h"
#include "./PrecMatrix.h"
#include "./NutMatrix.h"
#include "./PoleMatrix.h"
#include "./GHAMatrix.h"
#include "./AccelHarmonic.h"
#include "./G_AccelHarmonic.h"

//------------------------------------------------------------------------------
// VarEqn(double x,const Matrix& yPhi)
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

Matrix VarEqn(double x,const Matrix& yPhi);

#endif //VAREQN_H
