#include "../include/Matrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#define TOL_ 10e-14

Matrix::Matrix() : fil(0), col(0),matrix(nullptr){}

Matrix::Matrix(int fil, int col) : fil(fil), col(col)
{
    initMatrix();
}


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


Matrix::Matrix(const Matrix& m):Matrix()
{
    *this = m;
}


Matrix::~Matrix()
{
    for (int i = 0; i < fil; i++)
        delete[] matrix[i];
 
    delete[] matrix;
}


void Matrix::initMatrix()
{
    matrix = new double*[fil];
    for (int i = 0; i < fil; i++)
        matrix[i] = new double[col];
 
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            matrix[i][j] = 0.0;
}


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


Matrix Matrix::operator+(const Matrix& matrix2) const
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + matrix2.matrix[i][j];
 
    return result;
}


Matrix Matrix::operator+(double sumando) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] + sumando;

    return result;
}


Matrix operator+(const double& escalar,const Matrix& matrix)
{
    Matrix result(matrix.fil, matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] + escalar;

    return result;
}


Matrix Matrix::operator-(const Matrix& matrix2) const
{
    Matrix result(fil, col);
    
    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - matrix2.matrix[i][j];
 
    return result;
}


Matrix Matrix::operator-(double escalar) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] - escalar;

    return result;
}


Matrix operator-(const double& escalar,const Matrix& matrix)
{
    Matrix result(matrix.fil,matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] - escalar;

    return result;
}


Matrix Matrix::operator-()
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = -matrix[i][j];

    return result;
}

 
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


Matrix Matrix::operator*(double multiplicador) const
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] * multiplicador;

    return result;
}


Matrix operator*(const double &escalar, const Matrix &matrix) {
    Matrix result(matrix.fil, matrix.col);

    for (int i = 0; i < matrix.fil; i++)
        for (int j = 0; j < matrix.col; j++)
            result.matrix[i][j] = matrix.matrix[i][j] * escalar;

    return result;
}


Matrix Matrix::operator/(double divisor)
{
    Matrix result(fil, col);

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            result.matrix[i][j] = matrix[i][j] / divisor;

    return result;
}

 
double& Matrix::operator()(const int i, const int j) const
{
    return matrix[i-1][j-1];
}


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


Matrix Matrix::trans() const
{
    Matrix result(col,fil);

    for (int i = 0; i < col; i++)
        for (int j = 0; j < fil; j++)
            result.matrix[i][j] = matrix[j][i];

    return result;
}


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



int Matrix::getFil() const {
    return fil;
}


int Matrix::getCol() const {
    return col;
}


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


bool Matrix::equals(const Matrix &matrix2,const double TOL) const {
    if(fil != matrix2.fil || col != matrix2.col)
        return false;

    for (int i = 0; i < fil; i++)
        for (int j = 0; j < col; j++)
            if(fabs(matrix[i][j] - matrix2.matrix[i][j]) > TOL)
                return false;
    return true;
}


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


Matrix Matrix::eye(int n) {
    Matrix id(n,n);

    for (int i = 1; i <= n; i++)
        id(i,i) = 1;

    return id;
}


double Matrix::dot(const Matrix &matrix1, const Matrix &matrix2) {
    if(matrix1.col != matrix2.col)
        throw "Different dimensions";

    double sum = 0.0;

    for (int j = 1; j <= matrix1.col; j++){
        sum += matrix1(1,j)*matrix2(1,j);
    }

    return sum;
}


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
