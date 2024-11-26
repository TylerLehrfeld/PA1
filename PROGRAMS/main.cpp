#include <iostream>
#include <vector>

#include "DataReader.h"
#include "Matrix.h"
#include "Pivot.h"
#include "PointCloudTransform.h"
#include "Transform.h"
using std::cout;
using std::endl;

int main() {
    // Goal 1: develop a matrix module: run './test matrix' and see Matrix.h,
    // Matrix.cpp, and Matrix-test.h Goal 1: develop a transformation module:
    // run './test transform' and see Transform.h, Transform.cpp, and
    // Transform-test.cpp Goal 2: develop a 3D point set to point set
    // registration algorithm: run './test point-cloud' and see
    // PointCloudGenerator.h, PointCloudGenerator.cpp, and PointCloudTest.h Goal
    // 3: develop a pivot calibration method: run './test pivot' and see
    // Pivot.cpp, Pivot.h, and Pivot-test.h

    // for each debugging file set
    DataReader data = DataReader();
    for(int DEBUGGING_SET = 0; DEBUGGING_SET < data.numDataSets;
        DEBUGGING_SET++) {
        // read data
        data.reset();
        data.getCalBodyData(DEBUGGING_SET);
        data.getCalReadings(DEBUGGING_SET);
        data.getEMPivotReadings(DEBUGGING_SET);
        data.getLEDPivotReadings(DEBUGGING_SET);
        int N_Frames = data.calibrationDataFrames.size();
        int N_C = data.CalObjectEMPoints.size();
        string outputFileName = "OUTPUTFILE_" + std::to_string(DEBUGGING_SET);
        vector<Matrix> Cs;

        // Goal 4a: for each calibration data frame compute F_D
        PointCloudTransform T = PointCloudTransform();

        for(int TRANSFORM_NUM = 0;
            TRANSFORM_NUM < data.calibrationDataFrames.size();
            TRANSFORM_NUM++) {
            vector<Matrix> DCloud;
            vector<Matrix> ACloud;
            vector<Matrix> CCloud;
            std::tie(DCloud, ACloud, CCloud) =
                data.calibrationDataFrames[TRANSFORM_NUM];
            Transform F_D = T.compute(data.EMtrackerLEDPoints, DCloud);
            // Goal 4b: compute F_A
            Transform F_A = T.compute(data.CalObjectLEDPoints, ACloud);
            // Goal 4c: compute each C_i using F_D^-1*F_A*c_i
            for(int i = 0; i < data.CalObjectEMPoints.size(); i++) {
                Matrix C_i = F_D.inverse() * F_A * data.CalObjectEMPoints[i];
                Cs.push_back(C_i);
            }
            
        }
        // Goal 5: do a pivot calibration for EM tracking data
        // 5a
        Matrix G_0(3, 1, {0, 0, 0});
        for(int i = 0; i < data.N_G; i++) {
            G_0 = G_0 + data.EMPivotPointCloudFrames[0][i];
        }
        G_0 = 1 / (float)(data.N_G) * G_0;
        // 5b
        vector<Matrix> G_Hat;
        for(int i = 0; i < data.N_G; i++) {
            G_Hat.push_back(data.EMPivotPointCloudFrames[0][i] + -1 * G_0);
        }
        vector<Transform> F_Gs;
        for(int i = 0; i < data.N_GFrames; i++) {
            F_Gs.push_back(T.compute(G_Hat, data.EMPivotPointCloudFrames[i]));
        }
        // 5c
        Pivot p1 = Pivot(F_Gs);
        Matrix EMEst = p1.p_post;

        // Goal 6: Do a pivot calibration for LED data
        //repeat steps from 5
        Matrix H_0(3, 1, {0, 0, 0});
        vector<Matrix> D;
        vector<Matrix> H;
        std::tie(D, H) = data.LEDPivotPointCloudFrames[0];
        for(int i = 0; i < data.N_H; i++) {
            H_0 = H_0 + H[i];
        }

        H_0 = 1 / (float)(data.N_H) * H_0;
        vector<Matrix> H_Hat;
        for(int i = 0; i < data.N_H; i++) {
            H_Hat.push_back(H[i] + -1 * H_0);
        }
        // this time "first use F_D to transform the optical tracker"
        Transform F_D;
        vector<vector<Matrix>> F_DxH;
        for(int i = 0; i < data.N_HFrames; i++) {
            std::tie(D, H) = data.LEDPivotPointCloudFrames[i];
            //compute transform F_D
            F_D = T.compute(D, data.EMtrackerLEDPoints);
            vector<Matrix> F_DxH_i;
            for(int j = 0; j < H.size(); j++) {
                F_DxH_i.push_back(F_D * H[j]);
            }
            F_DxH.push_back(F_DxH_i);
        }
        vector<Transform> F_Hs;
        for(int i = 0; i < data.N_HFrames; i++) {
            F_Hs.push_back(T.compute(H_Hat, F_DxH[i]));
        }

        Pivot p2(F_Hs);
        Matrix OPTEst = p2.p_post;
        data.createOutPutFile(N_C, N_Frames, outputFileName, EMEst, OPTEst, Cs);
        if(DEBUGGING_SET < 8) {
            float error = data.averageError(EMEst, OPTEst, Cs, DEBUGGING_SET);
            // Goal 4d: output C_i expected
            cout << "Average error for C_i set " << DEBUGGING_SET << ": "
                 << error << endl;
        }
    }

    return 0;
}