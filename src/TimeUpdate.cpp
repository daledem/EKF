//$Header$
//
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/05/11
//
//------------------------------------------------------------------------------
#include "../include/TimeUpdate.h"

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix TimeUpdate(const Matrix &P,const Matrix &Phi, double Qdt)
//------------------------------------------------------------------------------
/**
 *   I don´t know what this does
 *
 * @param <P> Matrix
 * @param <Phi> Matrix
 * @param <Qdt> double
 *
 * @return Matrix
 */
//------------------------------------------------------------------------------
Matrix TimeUpdate(const Matrix &P,const Matrix &Phi, double Qdt) {
    return Phi*P*Phi.trans() + Qdt;
}
