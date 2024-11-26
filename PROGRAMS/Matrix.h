#include <vector>
#include <tuple>
#ifndef MATRIX
#define MATRIX

static const float MAX_ERROR = .00001;
static bool isCloseTo(float a, float b) {
    return (a - MAX_ERROR) < b && (a + MAX_ERROR) > b;
}

using std::vector;

class Matrix {
public:
    int ROWS;
    int COLUMNS;
    Matrix(int r, int c, vector<float> data);
    Matrix(vector<Matrix> columns);
    Matrix();
    Matrix inverse();
    ~Matrix();
    Matrix operator*(Matrix m2);
    Matrix operator+(Matrix m2);
    Matrix operator*(float scalar);
    Matrix adjunct();
    Matrix transpose();
    Matrix renormalize();
    Matrix cross(Matrix b);
    void switchRow(int row1, int row2);
    float divideRow(int row);
    void divideRow(int row, float divisor);
    float subtractRow(int row1, int row2, int column);
    void subtractRow(int row1, int row2, int column, float multiple);
    friend Matrix operator*(float scalar, Matrix m);
    bool operator==(Matrix b);
    float det();
    float magnitude();
    void print();
    
    std::tuple<vector<float>, vector<Matrix>> eigenValuesAndVectors();
    float trace();
    vector<float> matrixArray;
};

#endif