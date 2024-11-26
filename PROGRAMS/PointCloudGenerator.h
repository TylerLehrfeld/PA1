#include "Transform.h"
#include "Matrix.h"


#ifndef PCLOUDGEN
#define PCLOUDGEN
class PointCloudGenerator {
    public:
    PointCloudGenerator();
    vector<Matrix> PcloudA;
    Transform F_BA;
    vector<Matrix> PcloudB;
};


#endif