#include "Transform.h"

#include <stdexcept>

Transform::Transform(Transform AC, Transform CB) {
    R_AB = AC.R_AB * CB.R_AB;
    p_AB = AC.R_AB * CB.p_AB + AC.p_AB;
}

Matrix Transform::operator*(Matrix u_B) {
    return this->R_AB * u_B + p_AB;
}

Transform Transform::inverse() {
    Matrix inverse = (this->R_AB).inverse();
    return Transform((this->R_AB).inverse(), -1 * ((this->R_AB).inverse()*this->p_AB));
}

Transform::Transform() {
    R_AB = Matrix(3,3,{0,0,0,0,0,0,0,0,0});
    p_AB = Matrix(3,1,{0,0,0});
}

Transform::Transform(Matrix R, Matrix p) {
    if (R.COLUMNS != 3 || R.ROWS != 3) {
        throw std::invalid_argument(
            "rows and columns of rotation matrices must be 3.");
    }
    if (p.COLUMNS != 1 || p.ROWS != 3) {
        throw std::invalid_argument(
            "rows and columns of p vector must be 3 and 1 respectively.");
    }
    R_AB = R;
    p_AB = p;
}

Transform Transform::operator*(Transform T_BC) {
    return Transform(*this, T_BC);
}

void Transform::print() {
    R_AB.print();
    p_AB.print();
}
