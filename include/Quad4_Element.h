#ifndef QUAD4_ELEMENT_H
#define QUAD4_ELEMENT_H

#include <petscksp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <numeric>

#include "Element.h"

typedef Eigen::Array<PetscInt, -1, 1> ArrayXI;
typedef Eigen::Array<PetscInt, -1, -1> ArrayXXI;
typedef Eigen::Array<PetscInt, -1, -1, Eigen::RowMajor> ArrayXXIRM;
typedef Eigen::Array<PetscScalar, -1, -1, Eigen::RowMajor> ArrayXXSRM;
typedef Eigen::Array<PetscScalar, -1, 1> ArrayXSRM;


// Element class

class Quad4_Element: public Element
{
public:
    Quad4_Element();
public:
    void cal_ke();
};

#endif // ELEMENT_H
