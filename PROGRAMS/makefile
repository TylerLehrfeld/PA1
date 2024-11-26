# Compiler
CXX = g++

#paths for dependencies
PATHS = /usr/include/eigen3

# Flags
CXXFLAGS = -g -I $(PATHS)


# Program files
SOURCES = Matrix.cpp Transform.cpp PointCloudGenerator.cpp DataReader.cpp Pivot.cpp
HEADERS = Matrix.h Transform.h Matrix-test.h PointCloudTransform.h PointCloudTest.h PointCloudGenerator.h Transform-test.h Pivot.h Pivot-test.h
MAIN_OBJECTS = main.o $(SOURCES:.cpp=.o)
TEST_OBJECTS = test.o $(SOURCES:.cpp=.o)

# Main compilation
main: $(MAIN_OBJECTS)
	$(CXX) $(CXXFLAGS) $(MAIN_OBJECTS) -o main

# Test compilation
test: $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) $(TEST_OBJECTS) -o test

# Rule for compiling .cpp files into .o files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up object files and executables
clean:
	rm -f $(MAIN_OBJECTS) $(TEST_OBJECTS) main test Matrix-test.o

clean-nonexecutables:
	rm -f $(MAIN_OBJECTS) $(TEST_OBJECTS)  
