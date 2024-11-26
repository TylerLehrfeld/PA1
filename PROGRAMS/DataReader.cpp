#include "DataReader.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <cmath>
using std::cout, std::endl;
using std::unordered_map;

const vector<string> GroupLabels = {"a", "b", "c", "d", "e", "f",
                                    "g", "h", "i", "j", "k"};

vector<string> DataReader::getAllFileNames(string dir) {
    vector<string> filenames;
    for(auto const& file : std::filesystem::directory_iterator(dir)) {
        string filename = file.path().generic_string();
        filenames.push_back(filename);
    }
    return filenames;
}

void DataReader::createOutPutFile(int N_C, int N_F, string outputFileName,
                                  Matrix post1, Matrix post2,
                                  vector<Matrix> Cs) {
    cout << outputFileName << endl;
    std::ofstream outfile("../OUTPUT/"+outputFileName+".txt");
    if(!outfile.is_open()) {
        throw std::logic_error("The file wasn't opened");
    }
    outfile << N_C << ", " << N_F << ", " << outputFileName << "\n";
    outfile << post1.matrixArray[0] << ", " << post1.matrixArray[1] << ", " << post1.matrixArray[2] << "\n";
    outfile << post2.matrixArray[0] << ", " << post2.matrixArray[1] << ", " << post2.matrixArray[2] << "\n";
    for(int i = 0; i < Cs.size(); i++) {
        outfile << Cs[i].matrixArray[0] << ", " << Cs[i].matrixArray[1] << ", " << Cs[i].matrixArray[2] << "\n";    
    }
}


DataReader::DataReader() {
    vector<string> files = getAllFileNames("../data/PA 1 Student Data");
    organizeFiles(files);
}

void DataReader::organizeFiles(vector<string> files) {
    unordered_map<string, vector<string>> groupMap;
    while(files.size() != 0) {
        int last = files.size() - 1;
        for(int i = 0; i < GroupLabels.size(); i++) {
            if(static_cast<int>(files[last].find("-" + GroupLabels[i] + "-")) !=
               -1) {
                groupMap[GroupLabels[i]].push_back(files[last]);
                files.pop_back();
                i = GroupLabels.size();
            }
        }
    }
    for(int i = 0; i < GroupLabels.size(); i++) {
        groups.push_back(groupMap[GroupLabels[i]]);
    }
    numDataSets = groupMap.size();
}

void DataReader::reset() {
    EMtrackerLEDPoints.clear();
    CalObjectLEDPoints.clear();
    CalObjectEMPoints.clear();
    calibrationDataFrames.clear();
    EMPivotPointCloudFrames.clear();
    LEDPivotPointCloudFrames.clear();
}

float roundTo2Decimals(float num) {
    return std::ceil(num * 100.0) / 100.0;
}


float DataReader::averageError(Matrix p1, Matrix p2, vector<Matrix> cs, int groupIndex) {
    vector<string> files = groups[groupIndex];
    float error = 0;
    for(int i = 0; i < files.size(); i++) {
        if(files[i].find("output") != -1) {
            std::ifstream CalbodyFile(files[i]);
            string nextStr = "";
            CalbodyFile >> nextStr;
            int N_C = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            int N_Frames = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            string filename = nextStr;
            CalbodyFile >> nextStr;
            float p1xReal = std::stof(nextStr);
            float p1x = roundTo2Decimals(p1.matrixArray[0]);
            CalbodyFile >> nextStr;
            float p1yReal = std::stof(nextStr);
            float p1y = roundTo2Decimals(p1.matrixArray[1]);
            CalbodyFile >> nextStr;
            float p1zReal = std::stof(nextStr);
            float p1z = roundTo2Decimals(p1.matrixArray[2]);
            cout << "Magnitude of EMprobe error: " <<(Matrix(3,1,{p1xReal,p1yReal,p1zReal}) + -1* Matrix(3,1,{p1x,p1y,p1z})).magnitude() << endl;
            CalbodyFile >> nextStr;
            float p2xReal = std::stof(nextStr);
            float p2x = roundTo2Decimals(p2.matrixArray[0]);
            CalbodyFile >> nextStr;
            float p2yReal = std::stof(nextStr);
            float p2y = roundTo2Decimals(p2.matrixArray[1]);
            CalbodyFile >> nextStr;
            float p2zReal = std::stof(nextStr);
            float p2z = roundTo2Decimals(p2.matrixArray[2]);
            cout << "Magnitude of OpticalProbe error: " <<(Matrix(3,1,{p2xReal,p2yReal,p2zReal}) + -1* Matrix(3,1,{p2x,p2y,p2z})).magnitude() << endl;
            for(int i = 0; i < cs.size(); i++) {
                CalbodyFile >> nextStr;
                float d_Xi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float d_Yi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float d_Zi = std::stof(nextStr);
                Matrix c_i(3,1,{d_Xi,d_Yi,d_Zi});
                cs[i].matrixArray[0] = roundTo2Decimals(cs[i].matrixArray[0]);
                cs[i].matrixArray[1] = roundTo2Decimals(cs[i].matrixArray[1]);
                cs[i].matrixArray[2] = roundTo2Decimals(cs[i].matrixArray[2]);
                error += (cs[i] + -1 * c_i).magnitude();
            }
            
            
        }
    }
    return error/(float)cs.size();
}

void DataReader::getCalBodyData(int groupIndex) {
    vector<string> files = groups[groupIndex];
    for(int i = 0; i < files.size(); i++) {
        if(files[i].find("calbody") != -1) {
            std::ifstream CalbodyFile(files[i]);
            string nextStr = "";
            CalbodyFile >> nextStr;
            int N_D = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            int N_A = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            int N_C = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            string FileName = nextStr;
            for(int j = 0; j < N_D; j++) {
                CalbodyFile >> nextStr;
                float d_Xi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float d_Yi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float d_Zi = std::stof(nextStr);
                EMtrackerLEDPoints.push_back(Matrix(3, 1, {d_Xi, d_Yi, d_Zi}));
            }
            for(int j = 0; j < N_A; j++) {
                CalbodyFile >> nextStr;
                float a_Xi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float a_Yi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float a_Zi = std::stof(nextStr);
                CalObjectLEDPoints.push_back(Matrix(3, 1, {a_Xi, a_Yi, a_Zi}));
            }
            for(int j = 0; j < N_C; j++) {
                CalbodyFile >> nextStr;
                float c_Xi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float c_Yi = std::stof(nextStr);
                CalbodyFile >> nextStr;
                float c_Zi = std::stof(nextStr);
                CalObjectEMPoints.push_back(Matrix(3, 1, {c_Xi, c_Yi, c_Zi}));
            }
        }
    }
}

void DataReader::getCalReadings(int groupIndex) {
    vector<string> files = groups[groupIndex];
    for(int i = 0; i < files.size(); i++) {
        if(files[i].find("calreadings") != -1) {
            std::ifstream CalbodyFile(files[i]);
            string nextStr = "";
            CalbodyFile >> nextStr;
            int N_D = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            int N_A = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            int N_C = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            int N_Frames = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            string FileName = nextStr;
            for(int j = 0; j < N_Frames; j++) {
                vector<Matrix> DVectorCloud;
                for(int k = 0; k < N_D; k++) {
                    CalbodyFile >> nextStr;
                    float D_Xi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float D_Yi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float D_Zi = std::stof(nextStr);
                    DVectorCloud.push_back(Matrix(3, 1, {D_Xi, D_Yi, D_Zi}));
                }
                vector<Matrix> AVectorCloud;
                for(int k = 0; k < N_A; k++) {
                    CalbodyFile >> nextStr;
                    float A_Xi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float A_Yi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float A_Zi = std::stof(nextStr);
                    AVectorCloud.push_back(Matrix(3, 1, {A_Xi, A_Yi, A_Zi}));
                }
                vector<Matrix> CVectorCloud;
                for(int k = 0; k < N_C; k++) {
                    CalbodyFile >> nextStr;
                    float C_Xi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float C_Yi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float C_Zi = std::stof(nextStr);
                    CVectorCloud.push_back(Matrix(3, 1, {C_Xi, C_Yi, C_Zi}));
                }
                calibrationDataFrames.push_back(
                    std::make_tuple(DVectorCloud, AVectorCloud, CVectorCloud));
            }
            return;
        }
    }
}

void DataReader::getEMPivotReadings(int groupIndex) {
    vector<string> files = groups[groupIndex];
    for(int i = 0; i < files.size(); i++) {
        if(files[i].find("empivot") != -1) {
            std::ifstream CalbodyFile(files[i]);
            string nextStr = "";
            CalbodyFile >> nextStr;
            N_G = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            N_GFrames = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            string FileName = nextStr;
            for(int j = 0; j < N_GFrames; j++) {
                vector<Matrix> Gcloud;
                for(int k = 0; k < N_G; k++) {
                    CalbodyFile >> nextStr;
                    float G_Xi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float G_Yi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float G_Zi = std::stof(nextStr);
                    Gcloud.push_back(Matrix(3, 1, {G_Xi, G_Yi, G_Zi}));
                }
                EMPivotPointCloudFrames.push_back(Gcloud);
            }
            return;
        }
    }
}

void DataReader::getLEDPivotReadings(int groupIndex) {
    vector<string> files = groups[groupIndex];
    for(int i = 0; i < files.size(); i++) {
        if(files[i].find("optpivot") != -1) {
            std::ifstream CalbodyFile(files[i]);
            string nextStr = "";
            CalbodyFile >> nextStr;
            int N_D = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            N_H = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            N_HFrames = std::stoi(nextStr);
            CalbodyFile >> nextStr;
            string FileName = nextStr;
            for(int j = 0; j < N_HFrames; j++) {
                vector<Matrix> DCloud;
                for(int k = 0; k < N_D; k++) {
                    CalbodyFile >> nextStr;
                    float D_Xi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float D_Yi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float D_Zi = std::stof(nextStr);
                    DCloud.push_back(Matrix(3, 1, {D_Xi, D_Yi, D_Zi}));
                }
                vector<Matrix> HCloud;
                for(int k = 0; k < N_H; k++) {
                    CalbodyFile >> nextStr;
                    float H_Xi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float H_Yi = std::stof(nextStr);
                    CalbodyFile >> nextStr;
                    float H_Zi = std::stof(nextStr);
                    HCloud.push_back(Matrix(3, 1, {H_Xi, H_Yi, H_Zi}));
                }
                LEDPivotPointCloudFrames.push_back(std::tie(DCloud, HCloud));
            }
            return;
        }
        
    }
}
