//$Header$
//------------------------------------------------------------------------------
//                                  Matrix
//------------------------------------------------------------------------------
// EKF_GEOS3
//
// Author: David Ledesma
// Created: 2024/04/11
//
/**
 * This class provides a way to emulate matlab matrixes
 */
//------------------------------------------------------------------------------
#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#define TOL_ 10e-14

//---------------------------------
// public methods
//---------------------------------

//------------------------------------------------------------------------------
// Matrix()
//------------------------------------------------------------------------------
/**
 *   Construct Matrix object with 0 rows, 0 columns and matrix as nullptr
 */
//------------------------------------------------------------------------------
Matrix::Matrix() : fil(0), col(0),matrix(nullptr){}

//------------------------------------------------------------------------------
// Matrix(int fil, int col)
//------------------------------------------------------------------------------
/**
 *   Construct Matrix object with 0 in all its indexes
 *
 *   @param <fil> number of rows of the matrix
 *   @param <col> number of columns of the matrix
 */
//------------------------------------------------------------------------------
Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}

//------------------------------------------------------------------------------
// Matrix(int fil, int col, double v[], int n)
//------------------------------------------------------------------------------
/**
 *   Construct Matrix object that takes the values in a vector to
 *   populate its indexes
 *
 *   @param <fil> number of rows of the matrix
 *   @param <col> number of columns of the matrix
 *   @param <v> vector (double)
 *   @param <n> size of v
 */
//------------------------------------------------------------------------------
Matrix::Matrix(int fil, int col, double v[], int n): fil(fil), col(col)
{
    initMatrix();
 
    int k = 0;
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++){
            if (k < n)
                matrix[i][j] = v[k++];
            else
                matrix[i][j] = 0;
        }
}

//------------------------------------------------------------------------------
// Matrix(const Matrix& m)
//------------------------------------------------------------------------------
/**
 *  A copy constructor
 */
//------------------------------------------------------------------------------
Matrix::Matrix(const Matrix& m):Matrix()
{
    *this = m;
}

//------------------------------------------------------------------------------
// ~Matrix()
//------------------------------------------------------------------------------
/**
 *   Deletes the Matrix
 */
//------------------------------------------------------------------------------
Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}

//------------------------------------------------------------------------------
// Matrix& operator=(const Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 *   Assignment operator
 *
 *   @param <matrix2> Matrix object whose values to use to set "this" Matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix& Matrix::operator=(const Matrix& matrix2)
{
    if(matrix != nullptr) {
        for (int i = 0; i < fil; i++)
            delete[] matrix[i];
        delete[] matrix;
    }

    this->fil = matrix2.fil;
    this->col = matrix2.col;

    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            this->matrix[i][j] = matrix2.matrix[i][j];

    return *this;
}

//------------------------------------------------------------------------------
// Matrix operator+(const Matrix& matrix2) const
//------------------------------------------------------------------------------
/**
 *   Computes the sum between two Matrix objects
 *
 *   @param <matrix2> Matrix object to be add to "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator+(const Matrix& matrix2) const
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}

//------------------------------------------------------------------------------
// Matrix operator+(double sumando) const
//------------------------------------------------------------------------------
/**
 *   Computes the sum between a Matrix object and a double
 *
 *   @param <sumando> double to be add from "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator+(double sumando) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + sumando;

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator+(const double& escalar,const Matrix& matrix)
//------------------------------------------------------------------------------
/**
 *   Computes the sum between a Matrix object and a double
 *
 *   @param <escalar> double to be add to matrix
 *   @param <matrix> Matrix object to be add to escalar
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix operator+(const double& escalar,const Matrix& matrix)
{
    Matrix result(matrix.fil, matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] + escalar;

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator-(const Matrix& matrix2) const
//------------------------------------------------------------------------------
/**
 *   Computes the difference between two Matrix objects
 *
 *   @param <matrix2> Matrix object to be subtracted to "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator-(const Matrix& matrix2) const
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}

//------------------------------------------------------------------------------
// Matrix operator-(double escalar) const
//------------------------------------------------------------------------------
/**
 *   Computes the difference between a Matrix object and a double
 *
 *   @param <escalar> double to be subtracted from "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator-(double escalar) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - escalar;

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator-(const double& escalar,const Matrix& matrix)
//------------------------------------------------------------------------------
/**
 *   Computes the difference between a Matrix object and a double
 *
 *   @param <escalar> double to be subtracted to matrix
 *   @param <matrix> Matrix object to be subtracted to escalar
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix operator-(const double& escalar,const Matrix& matrix)
{
    Matrix result(matrix.fil,matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] - escalar;

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator-()
//------------------------------------------------------------------------------
/**
 *   Computes the opposite of a Matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator-()
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = -matrix[i][j];

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator*(const Matrix& matrix2) const
//------------------------------------------------------------------------------
/**
 *   Computes the product between two Matrix objects
 *
 *   @param <matrix2> Matrix object to be multiply by "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator*(const Matrix& matrix2) const
{
    if(col != matrix2.fil)
        throw "Wrong dimensions for matrix multiplication";

    Matrix result(fil, matrix2.col);
 
    for (int i = 0; i < this->fil ; i++){
        for (int j = 0; j < matrix2.col; j++){
            result.matrix[i][j] = 0;
            for (int k = 0; k < this->col; k++){
                result.matrix[i][j] = result.matrix[i][j] + this->matrix[i][k] * matrix2.matrix[k][j];
            }
        }
    }
 
    return result;
}

//------------------------------------------------------------------------------
// Matrix operator*(double multiplicador) const
//------------------------------------------------------------------------------
/**
 *   Computes the product between a Matrix object and a double
 *
 *   @param <multiplicador> doble to be multiply to "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator*(double multiplicador) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] * multiplicador;

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator*(const double &escalar, const Matrix &matrix)
//------------------------------------------------------------------------------
/**
 *   Computes the product between a Matrix object and a double
 *
 *   @param <escalar> doble to be multiply to "this" matrix object
 *   @param <matrix> matrix object to be multiply to "this" double
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix operator*(const double &escalar, const Matrix &matrix) {
    Matrix result(matrix.fil, matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] * escalar;

    return result;
}

//------------------------------------------------------------------------------
// Matrix operator/(double divisor)
//------------------------------------------------------------------------------
/**
 *   Computes the division between a Matrix object and a double
 *
 *   @param <divisor> double to divide "this" matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::operator/(double divisor)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] / divisor;

    return result;
}

//------------------------------------------------------------------------------
// double& operator()(const int i, const int j)
//------------------------------------------------------------------------------
/**
 *   Allow access to the values of the matrix emulating how matlab works
 *
 *   @param <i> index of the row
 *   @param <j> index of the column
 *
 *   @return double at (i,j)
 */
//------------------------------------------------------------------------------
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}

//------------------------------------------------------------------------------
// double determinant()
//------------------------------------------------------------------------------
/**
 *   Computes the determinant of "this" Matrix object
 *
 *   @return Determinant of "this" Matrix object
 *
 *   @note This code was taken from https://www.sanfoundry.com/cpp-program-find-inverse-matrix/
 */
//------------------------------------------------------------------------------
double Matrix::determinant() {
    double det = 0;
    double **pd = matrix;
    switch (fil)
    {
        case 1:
            det = pd[0][0];
            return det;
        case 2:
        {
            det = pd[0][0] * pd[1][1] - pd[0][1] * pd[1][0];
            return det;
        }
            break;
        case 3:
        {
            double a = pd[0][0];
            /***
            a b c
            d e f
            g h i
 
            a b c a b c
            d e f d e f
            g h i g h i
 
            // det (A) = aei + bfg + cdh - afh - bdi - ceg.
            ***/
            double b = pd[0][1];
            double c = pd[0][2];
            double d = pd[1][0];
            double e = pd[1][1];
            double f = pd[1][2];
            double g = pd[2][0];
            double h = pd[2][1];
            double i = pd[2][2];
            det = (a * e * i + b * f * g + c * d * h);
            det = det - a * f * h;
            det = det - b * d * i;
            det = det - c * e * g;
            return det;
        }
            break;
        case 4:
        {
            Matrix *temp[4];
            for (int i = 0; i < 4; i++)
                temp[i] = new Matrix(3, 3);
            for (int k = 0; k < 4; k++)
            {
                for (int i = 1; i < 4; i++)
                {
                    int j1 = 0;
                    for (int j = 0; j < 4; j++)
                    {
                        if (k == j)
                            continue;
                        temp[k]->matrix[i - 1][j1++]
                                        = this->matrix[i][j];
                    }
                }
            }
            det = this->matrix[0][0] * temp[0]->determinant()
                            - this->matrix[0][1] * temp[1]->determinant()
                            + this->matrix[0][2] * temp[2]->determinant()
                            - this->matrix[0][3] * temp[3]->determinant();
            return det;
        }
            break;
        case 5:
        {
            Matrix *temp[5];
            for (int i = 0; i < 5; i++)
                temp[i] = new Matrix(4, 4);
            for (int k = 0; k < 5; k++)
            {
                for (int i = 1; i < 5; i++)
                {
                    int j1 = 0;
                    for (int j = 0; j < 5; j++)
                    {
                        if (k == j)
                            continue;
                        temp[k]->matrix[i - 1][j1++]
                                        = this->matrix[i][j];
                    }
                }
            }
            det = this->matrix[0][0] * temp[0]->determinant()
                            - this->matrix[0][1] * temp[1]->determinant()
                            + this->matrix[0][2] * temp[2]->determinant()
                            - this->matrix[0][3] * temp[3]->determinant()
                            + this->matrix[0][4] * temp[4]->determinant();
            return det;
        }
        default:
        {
            int DIM = fil;
            Matrix **temp = new Matrix*[DIM];
            for (int i = 0; i < DIM; i++)
                temp[i] = new Matrix( DIM - 1, DIM - 1);
            for (int k = 0; k < DIM; k++)
            {
                for (int i = 1; i < DIM; i++)
                {
                    int j1 = 0;
                    for (int j = 0; j < DIM; j++)
                    {
                        if (k == j)
                            continue;
                        temp[k]->matrix[i - 1][j1++]
                                        = this->matrix[i][j];
                    }
                }
            }
            det = 0;
            for (int k = 0; k < DIM; k++)
            {
                if ((k % 2) == 0)
                    det = det + (this->matrix[0][k]
                                    * temp[k]->determinant());
                else
                    det = det - (this->matrix[0][k]
                                    * temp[k]->determinant());
            }
            for (int i = 0; i < DIM; i++)
                delete temp[i];
            delete[] temp;
            return det;
        }
            break;
    }
}

//------------------------------------------------------------------------------
// Matrix coFactor()
//------------------------------------------------------------------------------
/**
 *   Computes the coFactor of "this" Matrix object
 *
 *   @return coFator (as Matrix object) of "this" Matrix object
 *
 *   @note This code was taken from https://www.sanfoundry.com/cpp-program-find-inverse-matrix/
 */
//------------------------------------------------------------------------------
Matrix Matrix::coFactor() {
    Matrix cofactor(fil, col);
    if (fil != col)
        return cofactor;
    if (fil < 2)
        return Matrix::eye(fil);
    else if (fil == 2)
    {
        cofactor.matrix[0][0] = matrix[1][1];
        cofactor.matrix[0][1] = -matrix[1][0];
        cofactor.matrix[1][0] = -matrix[0][1];
        cofactor.matrix[1][1] = matrix[0][0];
        return cofactor;
    }
    else if (fil >= 3)
    {
        int DIM = fil;
        Matrix ***temp = new Matrix**[DIM];
        for (int i = 0; i < DIM; i++)
            temp[i] = new Matrix*[DIM];
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                temp[i][j] = new Matrix(DIM - 1, DIM - 1);
        for (int k1 = 0; k1 < DIM; k1++)
        {
            for (int k2 = 0; k2 < DIM; k2++)
            {
                int i1 = 0;
                for (int i = 0; i < DIM; i++)
                {
                    int j1 = 0;
                    for (int j = 0; j < DIM; j++)
                    {
                        if (k1 == i || k2 == j)
                            continue;
                        temp[k1][k2]->matrix[i1][j1++]
                                        = this->matrix[i][j];
                    }
                    if (k1 != i)
                        i1++;
                }
            }
        }
        bool flagPositive = true;
        for (int k1 = 0; k1 < DIM; k1++)
        {
            flagPositive = ((k1 % 2) == 0);
            for (int k2 = 0; k2 < DIM; k2++)
            {
                if (flagPositive == true)
                {
                    cofactor.matrix[k1][k2]
                                    = temp[k1][k2]->determinant();
                    flagPositive = false;
                }
                else
                {
                    cofactor.matrix[k1][k2]
                                    = -temp[k1][k2]->determinant();
                    flagPositive = true;
                }
            }
        }
        for (int i = 0; i < DIM; i++)
            for (int j = 0; j < DIM; j++)
                delete temp[i][j];
        for (int i = 0; i < DIM; i++)
            delete[] temp[i];
        delete[] temp;
    }
    return cofactor;
}

//------------------------------------------------------------------------------
// Matrix adjoint()
//------------------------------------------------------------------------------
/**
 *   Computes the adjoint of "this" Matrix object
 *
 *   @return adjoint (as Matrix object) of "this" Matrix object
 *
 *   @note This code was taken from https://www.sanfoundry.com/cpp-program-find-inverse-matrix/
 */
//------------------------------------------------------------------------------
Matrix Matrix::adjoint() {
    Matrix cofactor(fil, col);
    Matrix adj(fil, col);
    if (fil != col)
        return adj;
    cofactor = this->coFactor();
    // adjoint is transpose of a cofactor of a matrix
    for (int i = 0; i < fil; i++)
    {
        for (int j = 0; j < col; j++)
        {
            adj.matrix[j][i] = cofactor.matrix[i][j];
        }
    }
    return adj;
}

//------------------------------------------------------------------------------
// Matrix inverse()
//------------------------------------------------------------------------------
/**
 *   Computes the inverse of "this" Matrix object
 *
 *   @return inverse (as Matrix object) of "this" Matrix object
 *
 *   @note This code was taken from https://www.sanfoundry.com/cpp-program-find-inverse-matrix/
 */
//------------------------------------------------------------------------------
Matrix Matrix::inverse() {
    Matrix cofactor(fil, col);
    Matrix inv(fil, col);
    if (fil != col)
        return inv;
    // to find out Determinant
    double det = determinant();
    cofactor = this->coFactor();
    // inv = transpose of cofactor / Determinant
    for (int i = 0; i < fil; i++)
    {
        for (int j = 0; j < col; j++)
        {
            inv.matrix[j][i] = cofactor.matrix[i][j] / det;
        }
    }
    return inv;
}

//------------------------------------------------------------------------------
// Matrix trans() const
//------------------------------------------------------------------------------
/**
 *   Computes the transpose of "this" Matrix object
 *
 *   @return transpose (as Matrix object) of "this" Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::trans() const
{
    Matrix result(col,fil);

    for (int i = 0; i < col; i++)
        for (int j = 0; j < fil; j++)
            result.matrix[i][j] = matrix[j][i];

    return result;
}

//------------------------------------------------------------------------------
// Matrix append(Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 *   Appends the colums of two different Matrix objects with the same numeber of rows
 *
 *   @param <matrix2> Matrix object to be append to "this" Matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::append(Matrix& matrix2) {
    if(fil != matrix2.fil)
        throw "Different row dimensions";

    Matrix result(fil,col + matrix2.col);

    for (int i = 0; i < fil; i++) {
        int j;
        for (j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j];

        for (int k = 0; k < matrix2.col; j++,k++)
            result.matrix[i][j] = matrix2.matrix[i][k];
    }

    return result;
}

//------------------------------------------------------------------------------
// Matrix join(Matrix& matrix2)
//------------------------------------------------------------------------------
/**
 *   Joins the rows of two different Matrix objects with the same numeber of columns
 *
 *   @param <matrix2> Matrix object to be join to "this" Matrix object
 *
 *   @return Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::join(Matrix &matrix2) {
    if(col != matrix2.col)
        throw "Different col dimensions";

    Matrix result(fil + matrix2.fil,col);

    for (int i = 0; i < col; i++) {
        int j;
        for (j = 0; j < fil; j++)
            result.matrix[j][i] = matrix[j][i];

        for (int k = 0; k < matrix2.fil; j++,k++)
            result.matrix[j][i] = matrix2.matrix[k][i];
    }

    return result;
}

//------------------------------------------------------------------------------
// int getFil() const
//------------------------------------------------------------------------------
/**
 *   Gives the number of rows a Matrix object has
 *
 *   @return The number of rows
 */
//------------------------------------------------------------------------------
int Matrix::getFil() const {
    return fil;
}

//------------------------------------------------------------------------------
// int getCol() const
//------------------------------------------------------------------------------
/**
 *   Gives the number of columns a Matrix object has
 *
 *   @return The number of columns
 */
//------------------------------------------------------------------------------
int Matrix::getCol() const {
    return col;
}

//------------------------------------------------------------------------------
// Matrix getFilaByIndex(int fil,int inicio) const
//------------------------------------------------------------------------------
/**
 *   Gives the specify row from the specify start to the end
 *
 *   @param <fil> index of the row you want to obtain
 *   @param <inicio> index of the start (default is 0)
 *
 *   @return The row as a Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::getFilaByIndex(int fil,int inicio) const {
    if(inicio > col) {
        printf("Parameters outside of dimensions: %d",this->col);
        exit(EXIT_FAILURE);
    }

    Matrix result(1,col - inicio + 1);

    for (int i = inicio,j = 1; i <= col; i++,j++)
        result.matrix[0][j - 1] = matrix[fil - 1][i - 1];

    return result;
}

//------------------------------------------------------------------------------
// Matrix getFilaByIndex(int fil,int inicio, int fin) const
//------------------------------------------------------------------------------
/**
 *   Gives the specify row from the specify start to the specify end
 *
 *   @param <fil> index of the row you want to obtain
 *   @param <inicio> index of the start
 *   @param <fin> index of the end
 *
 *   @return The row as a Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::getFilaByIndex(int fil,int inicio, int fin) const {
    if(fin > col || inicio > col || inicio > fin) {
        printf("Parameters outside of dimensions: %d",this->col);
        exit(EXIT_FAILURE);
    }

    Matrix result(1,fin - inicio + 1);

    for (int i = inicio, j = 1; i <= fin; i++, j++)
        result.matrix[0][j - 1] = matrix[fil - 1][i - 1];

    return result;
}

//------------------------------------------------------------------------------
// Matrix getColumnaByIndex(int col,int inicio) const
//------------------------------------------------------------------------------
/**
 *   Gives the specify column from the specify start to the end
 *
 *   @param <col> index of the column you want to obtain
 *   @param <inicio> index of the start (default is 0)
 *
 *   @return The column as a Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::getColumnaByIndex(int col,int inicio) const {
    if(inicio > fil) {
        printf("Parameters outside of dimensions: %d",this->fil);
        exit(EXIT_FAILURE);
    }

    Matrix result(fil - inicio + 1, 1);

    for (int i = inicio; i <= fil; i++)
        result.matrix[i - 1][0] = matrix[i - 1][col - 1];

    return result;
}

//------------------------------------------------------------------------------
// Matrix getColumnaByIndex(int col,int inicio, int fin) const
//------------------------------------------------------------------------------
/**
 *   Gives the specify column from the specify start to the specify end
 *
 *   @param <col> index of the column you want to obtain
 *   @param <inicio> index of the start
 *   @param <fin> index of the end
 *
 *   @return The column as a Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::getColumnaByIndex(int col,int inicio,int fin) const {
    if(fin > fil || inicio > fil || inicio > fin) {
        printf("Parameters outside of dimensions: %d",this->fil);
        exit(EXIT_FAILURE);
    }

    Matrix result(fin - inicio + 1, 1);

    for (int i = inicio, j = 1; i <= fin; i++, j++)
        result.matrix[j - 1][0] = matrix[i - 1][col - 1];

    return result;
}

//------------------------------------------------------------------------------
// bool equals(const Matrix &matrix2,const double TOL) const
//------------------------------------------------------------------------------
/**
 *   Checks the similarity between two Matrix objects with a given tolerance
 *
 *   @param <matrix2> Matrix object to be compared to "this" Matrix object
 *   @param <TOL> The tolerance of difference that is allow between the two Matrix objects
 *
 *   @return true if they are equals with the given TOL, false otherwise
 */
//------------------------------------------------------------------------------
bool Matrix::equals(const Matrix &matrix2,const double TOL) const {
    if(fil != matrix2.fil || col != matrix2.col)
        return false;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            if(fabs(matrix[i][j] - matrix2.matrix[i][j]) > TOL)
                return false;
    return true;
}

//------------------------------------------------------------------------------
// void print()
//------------------------------------------------------------------------------
/**
 *   Prints the values in the Matrix object in the console
 */
//------------------------------------------------------------------------------
void Matrix::print()
{
    for (int i = 0; i < fil; i++){
        for (int j = 0; j < col; j++){
            std::cout << std::fixed << std::setprecision(14) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//------------------------------------------------------------------------------
// Matrix eye(int n)
//------------------------------------------------------------------------------
/**
 *   Creates an identity matrix of the given size
 *
 *   @param <n> size of the matrix
 *
 *   @return Identity matrix of size nxn
 */
//------------------------------------------------------------------------------
Matrix Matrix::eye(int n) {
    Matrix id(n,n);

    for (int i = 1; i <= n; i++)
        id(i,i) = 1;

    return id;
}

//------------------------------------------------------------------------------
// Matrix cross(const Matrix &matrix1, const Matrix &matrix2)
//------------------------------------------------------------------------------
/**
 *   Computes the cross product of two vectors
 *
 *   @param <matrix1> Matrix object with either 1 row or 1 column
 *   @param <matrix2> Matrix object with either 1 row or 1 column
 *
 *   @return The cross product as a Matrix object
 */
//------------------------------------------------------------------------------
Matrix Matrix::cross(const Matrix &matrix1, const Matrix &matrix2) {
    if(matrix1.fil == 1 && matrix1.col == 3 && matrix1.col == matrix2.col && matrix1.fil == matrix2.fil) {
        Matrix result(1,matrix1.col);

        result(1,1) = -matrix1(1,3)*matrix2(1,2) + matrix1(1,2)*matrix2(1,3);
        result(1,2) = matrix1(1,3)*matrix2(1,1) - matrix1(1,1)*matrix2(1,3);
        result(1,3) =  -matrix1(1,2)*matrix2(1,1) + matrix1(1,1)*matrix2(1,2);

        return result;
    }

    if(matrix1.col == 1 && matrix1.fil == 3 && matrix1.fil == matrix2.fil && matrix1.col == matrix2.col) {
        Matrix result(matrix1.fil,1);

        result(1,1) = -matrix1(3,1)*matrix2(2,1) + matrix1(2,1)*matrix2(3,1);
        result(2,1) = matrix1(3,1)*matrix2(1,1) - matrix1(1,1)*matrix2(3,1);
        result(3,1) =  -matrix1(2,1)*matrix2(1,1) + matrix1(1,1)*matrix2(2,1);

        return result;
    }

    printf("Wrong matrix dimensions for dot");
    exit(EXIT_FAILURE);
}

//------------------------------------------------------------------------------
// double dot(const Matrix &matrix1, const Matrix &matrix2)
//------------------------------------------------------------------------------
/**
 *   Computes the dot product of two vectors
 *
 *   @param <matrix1> Matrix object with either 1 row or 1 column
 *   @param <matrix2> Matrix object with either 1 row or 1 column
 *
 *   @return The dot product of the two vectors
 */
//------------------------------------------------------------------------------
double Matrix::dot(const Matrix &matrix1, const Matrix &matrix2) {
    double sum = 0.0;

    if(matrix1.fil == 1 && matrix1.col == matrix2.col && matrix1.fil == matrix2.fil) {
        for (int j = 1; j <= matrix1.col; j++){
            sum += matrix1(1,j)*matrix2(1,j);
        }

        return sum;
    }

    if(matrix1.col == 1 && matrix1.fil == matrix2.fil && matrix1.col == matrix2.col) {
        for (int j = 1; j <= matrix1.fil; j++){
            sum += matrix1(j,1)*matrix2(j,1);
        }

        return sum;
    }

    printf("Wrong matrix dimensions for dot");
    exit(EXIT_FAILURE);
}

//------------------------------------------------------------------------------
// double norm(const Matrix& matrix)
//------------------------------------------------------------------------------
/**
 *   Computes the norm of a vector
 *
 *   @param <matrix> Matrix object with either 1 row or 1 column
 *
 *   @return The norm of the vector
 */
//------------------------------------------------------------------------------
double Matrix::norm(const Matrix& matrix) {
    double sum = 0.0;
    if(matrix.fil == 1) {
        for (int j = 1; j <= matrix.col; j++){
            sum += matrix(1,j)*matrix(1,j);
        }

        return sqrt(sum);
    }
    if(matrix.col == 1) {
        for (int j = 1; j <= matrix.fil; j++){
            sum += matrix(j,1)*matrix(j,1);
        }

        return sqrt(sum);
    }

    printf("Wrong matrix dimensions for norm");
    exit(EXIT_FAILURE);

}

//---------------------------------
// private methods
//---------------------------------

//------------------------------------------------------------------------------
// void initMatrix()
//------------------------------------------------------------------------------
/**
 *   Allocates the space memory needed to save the values of the matrix
 */
//------------------------------------------------------------------------------
void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}