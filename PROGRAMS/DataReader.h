#include <vector>
#include <string>
#include "Matrix.h"
#include <tuple>
using std::string;
using std::vector;
using std::tuple;
#ifndef DATAREADER
#define DATAREADER

class DataReader {
    public:
    //initialize files
    DataReader();
    void organizeFiles(vector<string> files);
    void reset();
    int numDataSets;
    vector<vector<string>> groups;
    //file info getters
    //CALBODY.txt
    void getCalBodyData(int group);
    vector<Matrix> EMtrackerLEDPoints;
    vector<Matrix> CalObjectLEDPoints;
    vector<Matrix> CalObjectEMPoints;

    //CALREADINGS.txt
    void getCalReadings(int group);
    //list of data frames pcloud 1 is OpticalTrackerReadingsOfEMTrackerLEDPoints (D), 
    //pcloud 2 is OpticalTrackerReadingsOfCalObjLEDPoints (A)
    // and pcloud 3 is EMTrackerReadingsOfCalObjEMPoints (C);
    vector<tuple<vector<Matrix>,vector<Matrix>, vector<Matrix>>> calibrationDataFrames;
    
    //EMPIVOT.txt
    void getEMPivotReadings(int groupIndex);
    //list of point Clouds measured by EM tracker. Each point cloud is a list of G coordinates
    vector<vector<Matrix>> EMPivotPointCloudFrames;
    int N_G;
    int N_GFrames;

    //OPTPIVOT.txt
    void getLEDPivotReadings(int groupIndex);
    //list of point Frames containing two pointclouds: the first is the D cloud, the second is the H cloud
    vector<tuple<vector<Matrix>, vector<Matrix>>> LEDPivotPointCloudFrames;
    int N_H;
    int N_HFrames;
    void createOutPutFile(int N_C, int N_F, string outputFileName, Matrix post1, Matrix post2, vector<Matrix> Cs);
    float averageError(Matrix p1, Matrix p2, vector<Matrix> cs, int groupIndex);
    private:
    vector<string> getAllFileNames(string dir);

    };

#endif
