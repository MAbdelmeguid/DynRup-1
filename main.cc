#include <Eigen/Eigen>
#include <iostream>
#include "mesh.h"
int main()
{
	std::cout<<"hello world!"<<std::endl;
	Mesh Abaqus_io_mesh;

	string file_name = "test.inp";
	Abaqus_io_mesh.Abaqus_IO(file_name);



	std::cout<<Abaqus_io_mesh.Element.row(1)<<std::endl;


	return 1;


}
