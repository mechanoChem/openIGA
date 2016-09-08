#include "../../include/model/model.h"

/*
*root function to call all individual "residual" functions
*user may overwrite it
*/

template <class T, int dim>
void model<T, dim>::getResidual(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R, unsigned int currentIteration)
{
		//clean R check 
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
    if (abs(R[dof].val())>1.0e-13){
    	printf("**************Residual is contaminated**************. Value: %12.4e\n", R[dof].val());
			printf("\n"); exit(-1);
    }
	}
	iteration=currentIteration;
	
	residualForMechanics(cell, fe_values, ULocal, R);
	residualForNeummanBC(cell, fe_values, ULocal, R);
	  if (params->getBool("enforceWeakBC")) residualForHighOrderBC(cell,fe_values, ULocal, R);
}
template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;