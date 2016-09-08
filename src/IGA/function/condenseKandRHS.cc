/*
*IGA<dim>::condenseKandRHS
*Adjust right hand side and Jacobian (system_matrix) corresponding to dirrchelt BC
*/


#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::condenseKandRHS(){
  //Storing Residual value for output before condensing RHS
  // for (typename std::vector<controlPoint<dim> >::iterator i=mesh->controlPointVector.begin(); i<mesh->controlPointVector.end(); i++){
  // for (unsigned int a=0; a<dim; a++){
  // derivedValueR(i->id, a)=system_rhs(i->id*mesh->dofPerControlPoint + a);
  // }
  // }
  //Apply dirichlet map
  for (std::map<unsigned int, double>::iterator dof=dirichletMap.begin(); dof!=dirichletMap.end(); dof++){
    unsigned int i=dof->first;
    system_rhs(i)=0.0;
    for (unsigned int j=0; j<system_rhs.size(); j++){
      if ((i!=j) && (sparsity_pattern.nzMap.at(i).count(j)>0)) {system_matrix(i, j)=system_matrix(j, i)=0.0;}
    } 
  }
}

template class IGA<1>;
template class IGA<2>;
template class IGA<3>;
