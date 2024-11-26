#include <iostream>
#include <string>
#include <cmath>
#include "Matrix-test.h"
#include "PointCloudTest.h"
#include "Transform-test.h"
#include "Pivot-test.h"

using std::cout;
using std::endl;
using std::string;
int main(int argc, char *argv[]) {
    for(int i = 0; i < argc; i++) {
        string argument = argv[i];
        if(argument == "matrix") {
            testMatrixClass();
        }
        if(argument == "transform") {
            testTransformClass();
        }
        if(argument == "point-cloud") {
            testPointCloudClasses();
        }
        if(argument == "pivot") {
            testPivot();
        }
    }    
}