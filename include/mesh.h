#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

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

// Mesh class

class Mesh
{
public:
    Mesh(); 
    // ~Mesh() {std::cout << "Destructed\n";} // Desctuctor

    // MPI variables
    MPI_Comm comm;
    int myid;
    int nprocs;

    // Mesh information
    PetscInt nDims;
    ArrayXXSRM Node;
    ArrayXXIRM Element;
    ArrayXI gElem, gNode;
    ArrayXI elDist, ndDist;
    PetscInt nElem, nNode;
    PetscInt nLocElem, nLocNode;

    // FEM vectors
    Vec U, F;

    // Public functions
    PetscErrorCode Abaqus_IO(std::string &fname);
    PetscErrorCode Redistribute();

protected:
    // Auxiliary functions for redistributing mesh
    PetscErrorCode ReorderMETIS(PetscInt nparts=0, PetscInt ncommonNodes=0,
		                PetscScalar *tpwgts=NULL, PetscInt *elmwgt=NULL,
				PetscInt *opts=NULL);
    PetscErrorCode ReorderParMETIS(PetscInt nparts=0, PetscInt ncommonNodes=0,
		  PetscScalar *tpwgts=NULL, PetscScalar *ubvec=NULL,
		  PetscInt *opts=NULL, PetscInt ncon=1, PetscInt *elmwgt=NULL,
		  PetscInt wgtflag=0, PetscInt numflag=0 );
    PetscErrorCode ElemDist(Eigen::Array<PetscInt, -1, 1> &partition);
    PetscErrorCode NodeDist();
    PetscErrorCode Expand_Node();
    PetscErrorCode Initialize_Vectors();
    PetscErrorCode Localize();

};

#endif // MESH_H_INCLUDED
