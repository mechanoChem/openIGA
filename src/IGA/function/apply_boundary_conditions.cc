/*
*IGA<dim>::apply_boundary_conditions
*for user to apply dirichlet boundary condition
*/

#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::apply_boundary_conditions(){
  //Dirichlet map
  dirichletMap.clear(); 
  unsigned int  controlPointDOF=-dim;
  for (typename std::vector<controlPoint<dim> >::iterator controlpoint=mesh->controlPointVector.begin(); controlpoint<mesh->controlPointVector.end(); controlpoint++){
    std::vector<double> coords(controlpoint->coords);
    controlPointDOF+=dim;
    if( coords[1]==0.0){
      dirichletMap[controlPointDOF+0]=0.0;
      dirichletMap[controlPointDOF+1]=0.0;
      dirichletMap[controlPointDOF+2]=0.0;
    }

    //Apply values to solution vector
    for (std::map<unsigned int, double>::iterator dof=dirichletMap.begin(); dof!=dirichletMap.end(); dof++){
      U(dof->first)=dof->second;
    }
  }
	std::cout<<"finished BC"<<std::endl;
}

template class IGA<1>;
template class IGA<2>;
template class IGA<3>;