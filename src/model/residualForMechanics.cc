#include "../../include/model/model.h"
/*
*residual for regular stress and high order stress
*/
template <class T, int dim>
void model<T, dim>::residualForMechanics(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R){
  
  //evaluate Stress
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= fe_values.n_quadrature_points;
  dealii::Table<3,Sacado::Fad::DFad<double> > P (n_q_points, dim, dim);
  dealii::Table<4,Sacado::Fad::DFad<double> > Beta (n_q_points, dim, dim, dim);
  evaluateStress(fe_values, ULocal, P,Beta);
  
  //evaluate Residual
  for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
    R[dof] = 0;
    const unsigned int ck = fe_values.system_to_component_index(dof) - DOF;
    if (ck>=0 && ck<dim){
      // R = Grad(w)*P
      for (unsigned int q=0; q<n_q_points; ++q){
				for (unsigned int J = 0; J < dim; J++){
	  			R[dof] +=  fe_values.shape_grad(dof, q)[J]*P[q][ck][J]*fe_values.JxW(q);
 	 				for (unsigned int K = 0; K < dim; K++){
	    			R[dof] +=  fe_values.shape_grad_grad(dof, q)[J][K]*Beta[q][ck][J][K]*fe_values.JxW(q);
	  			}
				}
      }
    }
  }
}

template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;