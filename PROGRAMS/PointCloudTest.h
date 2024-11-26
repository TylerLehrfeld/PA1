#include <iostream>

#include "PointCloudGenerator.h"
#include "PointCloudTransform.h"
#include "helperFunctions.h"
using std::cout;
using std::endl;

void testPointCloudClasses() {
    cout << "testing generator" << endl;
    PointCloudTransform T = PointCloudTransform();
    for(int i = 0; i < 30; i++) {
        PointCloudGenerator p = PointCloudGenerator();
        Transform F_BA_Computed = T.compute(p.PcloudA, p.PcloudB);
        Transform IdentitiyTransform =
            p.F_BA * F_BA_Computed.inverse();
        assert(IdentitiyTransform.R_AB == I);
        assert(IdentitiyTransform.p_AB == origin);
        assert(p.F_BA * p.PcloudA[0] == F_BA_Computed * p.PcloudA[0]);
    }
    cout << "All point clouds are alligned: test point cloud succeeded" << endl;
}
