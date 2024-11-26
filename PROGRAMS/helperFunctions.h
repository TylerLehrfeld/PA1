#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "Matrix.h"
#include "Transform.h"
#include <cassert>

#ifndef HELPER
#define HELPER

static Matrix I(3, 3, {1, 0, 0, 0, 1, 0, 0, 0, 1});
static Matrix origin(3, 1, {0, 0, 0});

// generate a random float from 0-1
static float randomFloat() {
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

static Matrix generateRandomRotation() {

    float rotX = randomFloat() * 2 * M_PI;
    float rotY = randomFloat() * 2 * M_PI;
    float rotZ = randomFloat() * 2 * M_PI;

    // generate a rotation matrix using random angles
    Matrix Rx(
        3, 3,
        {1, 0, 0, 0, (float)cos(rotX), -sinf(rotX), 0, sinf(rotX), cosf(rotX)});
    Matrix Ry(3, 3,
              {cosf(rotY), 0, sinf(rotY), 0, 1, 0, -sinf(rotY), 0, cosf(rotY)});
    Matrix Rz(3, 3,
              {cosf(rotZ), -sinf(rotZ), 0, sinf(rotZ), cosf(rotZ), 0, 0, 0, 1});
    Matrix R = Rx * Ry * Rz;
    return R;
}

static Matrix generateRandomPoint() {
    // generate a random vector with values from -10 to +10
    Matrix p(3, 1,
             {randomFloat() * 20 - 10, randomFloat() * 20 - 10,
              randomFloat() * 20 - 10});
    return p;
}

static Transform generateRandomTransform() {
    // create the transform
    return Transform(generateRandomRotation(), generateRandomPoint());
}

static vector<Transform> generatePivotFrames(Matrix p_tip, Matrix p_post) {
    vector<Transform> pivots;
    // generate 15 frames of pivots
    for(int i = 0; i < 15; i++) {
        // first come up with a rotation, then work to get the displacement
        // based on the tip and post
        Matrix R_i = generateRandomRotation();
        // given p_post = R_i*p_tip+p_i, p_i = p_post - R_i*p_tip
        Matrix p_i = p_post + -1 * R_i * p_tip;
        pivots.push_back(Transform(R_i, p_i));
        assert(p_post == pivots[i] * p_tip);
    }
    return pivots;
}

#endif