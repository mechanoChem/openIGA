#include "../../include/model/model.h"
/*
*residual for neumman boundary condition
*user should overwrite it
*/
template <class T, int dim>
void model<T, dim>::residualForNeummanBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{
	load= params->getDouble("load"); 
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
//example for Newmman boundary condition; no need for defect paper.
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
    if (cell.boundaryFlags[faceID]>0){      
      //loop over boundary faces
      for (unsigned int dof=0; dof<dofs_per_cell; ++dof) {
				const unsigned int ck = fe_values.system_to_component_index(dof) - DOF;
				for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  			if (std::strcmp(bcType,"bending")==0){
	    			//bending along dim-1 direction
	    			if ((cell.boundaryFlags[faceID]==2*dim) and (ck==dim-2)){
	      			R[dof] += -fe_values.shape_value_face(dof, q, faceID)*load*fe_values.JxW_face(q, faceID);
	    			}
	  			}
	  		}
			} 
    } 
  }
}

template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;