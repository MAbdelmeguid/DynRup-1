#include <Eigen/Eigen>
#include <iostream>
#include "mesh.h"
#include "Element.h"
#include "Quad4_Element.h"

using namespace std;
static char help[] = "Is something breaking? Is it also breaking over time? Then this is the software for you.\n\n";

int main(int argc, char **args)
{
	PetscInt err = PetscInitialize(&argc, &args, (char*)0, help);
	std::cout<<"hello world!"<<std::endl;
	Mesh Abaqus_io_mesh;

	string file_name = "test.inp";
	Abaqus_io_mesh.Abaqus_IO(file_name);
	err = Abaqus_io_mesh.Redistribute(); CHKERRQ(err);

	err = PetscFinalize();

	// Temp test for Quad4 routine being called
	Quad4_Element Q4;
	Q4.cal_ke();

	return 0;


}
