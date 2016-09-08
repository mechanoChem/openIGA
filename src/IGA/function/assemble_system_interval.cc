/*
*IGA<dim>::assemble_system_interval
*for user to assemble system, usually no need to change except the model to call getesidual
*
*/

#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end){
  //element loop
  IGAValues<dim> fe_values_base(mesh, dim, 2);
  for (typename std::vector<knotSpan<dim> >::iterator cell=begin; cell<end; cell++){
    fe_values_base.reinit(*cell);
    IGAValues<dim>& fe_values=fe_values_base;
    //IGAValues<dim>* fe_values=cellValues[cell->id];
    unsigned int n_q_points= fe_values.n_quadrature_points;
    unsigned int dofs_per_cell=fe_values.dofs_per_cell;
    denseMatrix local_matrix(dofs_per_cell, dofs_per_cell);
    denseVector local_rhs(dofs_per_cell);
    //AD variables
    dealii::Table<1, Sacado::Fad::DFad<double> > ULocal(dofs_per_cell); dealii::Table<1, double > ULocalConv(dofs_per_cell);
    dealii::Table<1, double> ULocalTemp(dofs_per_cell); 
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      ULocal[i]=U(cell->local_dof_indices[i]);
      ULocalTemp[i]=U(cell->local_dof_indices[i]);
      ULocal[i].diff (i, dofs_per_cell);
      ULocalConv[i]= Un(cell->local_dof_indices[i]);
    }
	
    dealii::Table<1, Sacado::Fad::DFad<double> > R(dofs_per_cell); //R.fill(0.0);
    
    for(unsigned int i=0; i<dofs_per_cell; ++i) R[i]=0.0; 
		
		modelexample->getResidual(*cell, fe_values, ULocal, R, currentIteration);
			
    //Residual(R) and Jacobian(R')
    for (unsigned int i=0; i<dofs_per_cell; ++i) {
      for (unsigned int j=0; j<dofs_per_cell; ++j){
				// R' by AD
				local_matrix(i,j)= R[i].fastAccessDx(j);
      }
      local_rhs(i) = -R[i].val();
    }
	
    //Global Assembly
    assembler_lock.acquire();
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      for (unsigned int j=0; j<dofs_per_cell; ++j){
				system_matrix(cell->local_dof_indices[i], cell->local_dof_indices[j])+=local_matrix(i,j);
      }
      system_rhs(cell->local_dof_indices[i]) += local_rhs(i);
    }
    assembler_lock.release();
  }
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;

