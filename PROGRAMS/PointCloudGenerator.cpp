#include "PointCloudGenerator.h"

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include "helperFunctions.h"
using std::cout;
using std::endl;

PointCloudGenerator::PointCloudGenerator() {
    F_BA = generateRandomTransform();
    // create a list of points in 3D space (3x1 matrices);
    // then calculate its corresponding point in 3D space after transformation.
    int numPoints = (int)(randomFloat() * 20 + 10);
    for(int i = 0; i < numPoints; i++) {
        PcloudA.push_back(
            Matrix(3, 1,
                   {randomFloat() * 20 - 10, randomFloat() * 20 - 10,
                    randomFloat() * 20 - 10}));
        PcloudB.push_back(F_BA * PcloudA[i]);
    }
}
