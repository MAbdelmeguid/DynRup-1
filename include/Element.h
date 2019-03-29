#ifndef ELEMENT_H
#define ELEMENT_H

#include <petscksp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Eigen>
#include <numeric>

typedef Eigen::Array<PetscInt, -1, 1> ArrayXI;
typedef Eigen::Array<PetscInt, -1, -1> ArrayXXI;
typedef Eigen::Array<PetscInt, -1, -1, Eigen::RowMajor> ArrayXXIRM;
typedef Eigen::Array<PetscScalar, -1, -1, Eigen::RowMajor> ArrayXXSRM;
typedef Eigen::Array<PetscScalar, -1, 1> ArrayXSRM;


// Element class

class Element
{
public:
    Element();

protected:
    virtual void cal_ke()=0;
    ArrayXXSRM ke;
};

#endif // ELEMENT_H
