#include "Matrix.h"
#include <Eigen/Dense>
#include <iostream>

#include "helperFunctions.h"
using std::cout;
using std::endl;
using std::vector;

Matrix::Matrix(int rows, int columns, vector<float> data) {
    ROWS = rows;
    COLUMNS = columns;
    matrixArray = data;
}

Matrix::Matrix(vector<Matrix> columns) {
    ROWS = columns[0].ROWS;
    COLUMNS = columns.size();
    matrixArray = vector<float>(ROWS * COLUMNS, 0);
    for(int i = 0; i < ROWS; i++) {
        for(int j = 0; j < COLUMNS; j++) {
            matrixArray[i * COLUMNS + j] = columns[j].matrixArray[i];
        }
    }
}

Matrix::Matrix() {
    ROWS = 3;
    COLUMNS = 3;
    matrixArray = {0, 0, 0, 0, 0, 0, 0, 0, 0};
}
Matrix Matrix::inverse() {
    if(ROWS != COLUMNS) {
        throw std::invalid_argument(
            "Must be a square matrix to take the determinant.");
    }
    if(ROWS == 2 || ROWS == 3) {
        float det = this->det();

        if(det == 0) {
            throw std::invalid_argument(
                "This matrix is not inversible: determinant is 0");
        }

        if(ROWS == 2) {
            return (1 / det) * Matrix(2, 2,
                                      {matrixArray[3], -matrixArray[1],
                                       -matrixArray[2], matrixArray[0]});
        }

        if(ROWS == 3) {
            return this->adjunct() * (1 / det);
        }
    }
    if(ROWS == 6) {
        vector<float> CopiedArray;
        for(int i = 0; i < 36; i++) {
            CopiedArray.push_back(matrixArray[i]);
        }
        Matrix inv(6, 6, CopiedArray);
        Matrix Identity(6, 6,
                        {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                         0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1});
        for(int i = 0; i < 6; i++) {
            // search for value in column I with non 0 value
            int row = i;
            while(isCloseTo(inv.matrixArray[row * 6 + i], 0)) {
                row++;
                if(row > 5) {
                    throw std::invalid_argument(
                        "this matrix is not inversible");
                }
            }
            // switch given row with row i in this matrix and Identity
            inv.switchRow(row, i);
            Identity.switchRow(row, i);
            // make it so that row i has a 1 in column i
            float divisor = inv.divideRow(i);
            Identity.divideRow(i, divisor);
            
            // subtract row i from each other row
            for(int j = 0; j < 6; j++) {
                if(i != j) {
                    float mult = inv.subtractRow(i, j, i);
                    Identity.subtractRow(i, j, i, mult);
                }
            }
        }
        assert(inv.matrixArray ==
               vector<float>({1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                              0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                              0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}));
        assert(
            *this * Identity ==
            Matrix(6, 6, vector<float>({1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                        0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1})));
        return Identity;
    }
    throw std::domain_error(
        "inverse is only implemented for 2x2, 3x3, and 6x6 matrices.");
}

Matrix::~Matrix() {
}

float Matrix::det() {
    if(COLUMNS != ROWS) {
        throw std::invalid_argument(
            "Must be a square matrix to take the determinant.");
    }
    if(COLUMNS == 2) {
        return matrixArray[0] * matrixArray[3] -
               matrixArray[1] * matrixArray[2];
    } else if(COLUMNS == 3) {
        return matrixArray[0] * (matrixArray[4] * matrixArray[8] -
                                 matrixArray[5] * matrixArray[7]) -
               matrixArray[1] * (matrixArray[3] * matrixArray[8] -
                                 matrixArray[5] * matrixArray[6]) +
               matrixArray[2] * (matrixArray[3] * matrixArray[7] -
                                 matrixArray[4] * matrixArray[6]);
    } else {
        throw std::domain_error("Determinants of matrices other than 2x2 and "
                                "3x3 are not implemented.");
    }
}

float Matrix::magnitude() {
    if(COLUMNS != 1) {
        throw std::domain_error(
            "Magnitude not applicable to non column vectors");
    }
    float sum = 0;
    for(int i = 0; i < ROWS; i++) {
        sum += matrixArray[i] * matrixArray[i];
    }
    return std::sqrt(sum);
}

void Matrix::print() {
    for(int i = 0; i < ROWS; i++) {
        for(int j = 0; j < COLUMNS; j++) {
            cout << matrixArray[i * COLUMNS + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
std::tuple<vector<float>, vector<Matrix>> Matrix::eigenValuesAndVectors() {
    if(COLUMNS != 3 || ROWS != 3) {
        throw std::domain_error(
            "only implemented eigenvalues for 3x3 matrices");
    }
    Eigen::Matrix3d M;
    M << matrixArray[0], matrixArray[1], matrixArray[2], matrixArray[3],
        matrixArray[4], matrixArray[5], matrixArray[6], matrixArray[7],
        matrixArray[8];
    Eigen::EigenSolver<Eigen::Matrix3d> solver(M);

    // Get the eigenvalues
    Eigen::Vector3d eigenvalues = solver.eigenvalues().real();
    Eigen::Matrix3d eigenvectors = solver.eigenvectors().real();
    vector<Matrix> eigenvectorList;
    for(int i = 0; i < eigenvectors.cols(); i++) {
        eigenvectorList.push_back(Matrix(
            3, 1,
            vector<float>(
                eigenvectors.col(i).data(),
                eigenvectors.col(i).data() + eigenvectors.col(i).size())));
    }
    return std::make_tuple(
        vector<float>(eigenvalues.data(),
                      eigenvalues.data() + eigenvalues.size()),
        eigenvectorList);
}

float Matrix::trace() {
    if(COLUMNS != ROWS) {
        throw std::invalid_argument("rows must equal columns for trace");
    }
    float tr = 0;
    for(int i = 0; i < ROWS; i++) {
        tr += matrixArray[i * ROWS + i];
    }
    return tr;
}

Matrix Matrix::operator*(Matrix m2) {
    if(this->COLUMNS != m2.ROWS) {
        throw std::invalid_argument(
            "first array columns must equal second array rows");
    }

    vector<float> newData;
    // for each element in the new matrix
    for(int i = 0; i < this->ROWS; i++) {
        for(int j = 0; j < m2.COLUMNS; j++) {
            // calculate the new value of the matrix
            float dotProd = 0;
            for(int k = 0; k < this->COLUMNS; k++) {
                dotProd += this->matrixArray[this->COLUMNS * i + k] *
                           m2.matrixArray[k * m2.COLUMNS + j];
            }
            newData.push_back(dotProd);
        }
    }
    Matrix newMatrix(this->ROWS, m2.COLUMNS, newData);
    return newMatrix;
}

Matrix Matrix::operator+(Matrix m2) {
    if(this->COLUMNS != m2.COLUMNS || this->ROWS != m2.ROWS) {
        throw std::invalid_argument("The size of the matrices must match in "
                                    "order to add them together.");
    }
    vector<float> newMatrixValues;
    for(int i = 0; i < ROWS * COLUMNS; i++) {
        newMatrixValues.push_back(this->matrixArray[i] + m2.matrixArray[i]);
    }
    return Matrix(ROWS, COLUMNS, newMatrixValues);
}

Matrix Matrix::operator*(float scalar) {
    vector<float> newMatrixValues;
    for(int i = 0; i < ROWS * COLUMNS; i++) {
        newMatrixValues.push_back(this->matrixArray[i] * scalar);
    }
    return Matrix(this->ROWS, this->COLUMNS, newMatrixValues);
}

Matrix Matrix::adjunct() {
    if(ROWS != 3 || COLUMNS != 3) {
        throw std::invalid_argument(
            "The adjoint can only be taken for 3x3 matrices");
    }
    vector<float> adjointValues;
    adjointValues.push_back(
        Matrix(2, 2,
               {matrixArray[4], matrixArray[5], matrixArray[7], matrixArray[8]})
            .det());
    adjointValues.push_back(-Matrix(2, 2,
                                    {matrixArray[1], matrixArray[2],
                                     matrixArray[7], matrixArray[8]})
                                 .det());
    adjointValues.push_back(
        Matrix(2, 2,
               {matrixArray[1], matrixArray[2], matrixArray[4], matrixArray[5]})
            .det());
    adjointValues.push_back(-Matrix(2, 2,
                                    {matrixArray[3], matrixArray[5],
                                     matrixArray[6], matrixArray[8]})
                                 .det());
    adjointValues.push_back(
        Matrix(2, 2,
               {matrixArray[0], matrixArray[2], matrixArray[6], matrixArray[8]})
            .det());
    adjointValues.push_back(-Matrix(2, 2,
                                    {matrixArray[0], matrixArray[2],
                                     matrixArray[3], matrixArray[5]})
                                 .det());
    adjointValues.push_back(
        Matrix(2, 2,
               {matrixArray[3], matrixArray[4], matrixArray[6], matrixArray[7]})
            .det());
    adjointValues.push_back(-Matrix(2, 2,
                                    {matrixArray[0], matrixArray[1],
                                     matrixArray[6], matrixArray[7]})
                                 .det());
    adjointValues.push_back(
        Matrix(2, 2,
               {matrixArray[0], matrixArray[1], matrixArray[3], matrixArray[4]})
            .det());
    return Matrix(3, 3, adjointValues);
}

Matrix Matrix::transpose() {
    vector<float> transposedValues;
    for(int i = 0; i < COLUMNS; i++) {
        for(int j = 0; j < ROWS; j++) {
            transposedValues.push_back(matrixArray[j * COLUMNS + i]);
        }
    }
    return Matrix(COLUMNS, ROWS, transposedValues);
}

Matrix Matrix::renormalize() {
    if(ROWS != 3 || COLUMNS != 3) {
        throw std::invalid_argument("Cannot cross non 3x1 matrices");
    }
    Matrix r_y(3, 1, {matrixArray[1], matrixArray[4], matrixArray[7]});
    Matrix r_z(3, 1, {matrixArray[2], matrixArray[5], matrixArray[8]});
    Matrix a = r_y.cross(r_z);
    Matrix b = r_z.cross(a);
    a = 1 / ((a.transpose() * a).matrixArray[0]) * a;
    b = 1 / ((b.transpose() * b).matrixArray[0]) * b;
    r_z = 1 / ((r_z.transpose() * r_z).matrixArray[0]) * r_z;
    return Matrix(3, 3,
                  {a.matrixArray[0], b.matrixArray[0], r_z.matrixArray[0],
                   a.matrixArray[1], b.matrixArray[1], r_z.matrixArray[1],
                   a.matrixArray[2], b.matrixArray[2], r_z.matrixArray[2]});
}

Matrix Matrix::cross(Matrix b) {
    if(ROWS != 3 || b.ROWS != 3 || COLUMNS != 1 || b.COLUMNS != 1) {
        throw std::invalid_argument("Cannot cross non 3x1 matrices");
    }
    return Matrix(
        3, 1,
        {matrixArray[1] * b.matrixArray[2] - matrixArray[2] * b.matrixArray[1],
         matrixArray[2] * b.matrixArray[0] - matrixArray[0] * b.matrixArray[2],
         matrixArray[0] * b.matrixArray[1] -
             matrixArray[1] * b.matrixArray[0]});
}

Matrix operator*(float scalar, Matrix m) {
    return m * scalar;
}

bool Matrix::operator==(Matrix b) {
    if(COLUMNS != b.COLUMNS || ROWS != b.ROWS) {
        return false;
    }
    if(matrixArray.size() != b.matrixArray.size()) {
        return false;
    }
    for(size_t i = 0; i < matrixArray.size(); i++) {
        if(!isCloseTo(matrixArray[i], b.matrixArray[i])) {
            return false;
        }
    }

    return true;
}

void Matrix::switchRow(int row1, int row2) {
    for(int i = 0; i < COLUMNS; i++) {
        float temp = matrixArray[row2 * COLUMNS + i];
        matrixArray[row2 * COLUMNS + i] = matrixArray[row1 * COLUMNS + i];
        matrixArray[row1 * COLUMNS + i] = temp;
    }
}

float Matrix::divideRow(int row) {
    float divisor = matrixArray[row * COLUMNS + row];
    for(int i = 0; i < COLUMNS; i++) {
        matrixArray[row * COLUMNS + i] /= divisor;
    }
    return divisor;
}

void Matrix::divideRow(int row, float divisor) {
    for(int i = 0; i < COLUMNS; i++) {
        matrixArray[row * COLUMNS + i] /= divisor;
    }
}

float Matrix::subtractRow(int row1, int row2, int column) {
    float multiple = matrixArray[row2 * COLUMNS + column];
    for(int i = 0; i < COLUMNS; i++) {
        matrixArray[row2 * COLUMNS + i] -=
            matrixArray[row1 * COLUMNS + i] * multiple;
    }
    return multiple;
}

void Matrix::subtractRow(int row1, int row2, int column, float multiple) {
    for(int i = 0; i < COLUMNS; i++) {
        matrixArray[row2 * COLUMNS + i] -= multiple * matrixArray[row1 * COLUMNS + i];
    }
}
