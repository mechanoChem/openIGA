/*
*IGA<dim>::output
*Fill output vector and output vtk file
*/

#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::output (unsigned int _cycle){
	solutionClass<dim> displacement(*mesh, NODAL, VECTOR, std::string("u"));
	//fill output vector
  for (typename std::vector<controlPoint<dim> >::iterator i=mesh->controlPointVector.begin(); i<mesh->controlPointVector.end(); i++){
    for (unsigned int a=0; a<dim; a++){
			displacement(i->id, a)=Un(i->id*mesh->dofPerControlPoint + a);
    }
  }
	outputVariables.push_back(&displacement);
  char fileName[200];
  std::sprintf (fileName, "output0");
  std::cout<<"begin write mesh"<<std::endl;
  std::vector<int> outputGridSize(dim,10);
   writeMesh<dim>(fileName, _cycle, mesh, outputGridSize, outputVariables);
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;