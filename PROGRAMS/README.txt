DataReader.cpp: Reads data from data/PA 1 Student Data and writes Output files.
DataReader.h: Defines methods for DataReader
Dockerfile: defines environment if you cannot intall dependencies
helperFunctions.h: a set of helper functions and useful constants used in multiple files
main.cpp: calculates output files and puts them in ./OUTPUT dir with format OUTPUT_i
makefile: can either 'make' for main, or 'make test' for tests.
Matrix-test.h: a file of matrix tests for matrix functions
Matrix.cpp: a file implementing matrix functions
Matrix.h: a file defining the Matrix class.
Pivot-test.h: defines testPivot function which random pivots and then checks if the algorithm can make sense of them
Pivot.cpp: implementation for Pivot class
Pivot.h: define the Pivot class
PointCloudGenerator.cpp: generates a random point cloud for testing
PointCloudGenerator.h: defines the pointCloudGenerator class
PointCloudTest.h: definines test for point cloud registration
PointCloudTransform.h: implements computing a point cloud registration with SVD
test.cpp: runs all of the tests if given arguments: matrix transform pivot point-cloud
Transform-test.h: defines tests for Frame transformations
Transform.cpp: implements Transform functionality
Transform.h: definies the transform class



main: executable
test: executable

running steps:
1. go in PROGRAMS dir
2. run 'docker build -t matrix:first .'
3. cd ..
4. run 'docker run -it --mount "type=bind,source=$(pwd)/,target=/Homework1" matrix:first'
5. within the docker run 'cd Homework1/PROGRAMS' the './main' or './test matrix pivot point-cloud transform'


If you don't want to use a docker container, you may have to download Eigen yourself and configure the path.