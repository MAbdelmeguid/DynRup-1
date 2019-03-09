// C++ Standard Lib
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <numeric>
// Eigen Lib


// Mesh class

using namespace Eigen;
using namespace std;
class Mesh
{
public:
    // Mesh() {std::cout << "Constructed\n";} // Constructor
    // ~Mesh() {std::cout << "Destructed\n";} // Desctuctor

    ArrayXXd Node;
    ArrayXXi Element;

    void Abaqus_IO(string &fname);


};
