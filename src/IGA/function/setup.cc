/*
*IGA<dim>::setup
*Initialize vector and matrix
*/

#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::setup (){	
  //initialize global data structures
  sparsity_pattern.init(mesh);
  system_matrix.reinit(sparsity_pattern);
  system_rhs.reinit(sparsity_pattern); 
	U.reinit(sparsity_pattern); 
	Un.reinit(sparsity_pattern); 
	dU.reinit(sparsity_pattern); 
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;