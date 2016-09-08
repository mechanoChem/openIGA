#include "../../include/model/model.h"
/*
*Evaluate functions using quadrature points (first three dof e.g. displacement) interpolation
*Include scalar/vector functions and its derivatives 
*/
template <class T, int dim>
void model<T, dim>::evaluateScalarFunction(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<1, T>& U,  int faceID){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    U[q]=0.0; //U
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      if (fe_values.system_to_component_index(k)==DOF){
				if (faceID>-1){
	  			U[q]+=ULocal[k]*fe_values.shape_value_face(k, q, faceID); //U
				}
				else{
	  			U[q]+=ULocal[k]*fe_values.shape_value(k, q); 
				}
      }
    }
  }
}

template <class T, int dim>
void model<T, dim>::evaluateScalarFunctionGradient(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<2, T>& gradU, deformationMap<T, dim>& defMap, bool gradientInCurrentConfiguration, int faceID){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  dealii::Table<1, T> refGradU(dim);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    //refGradU.fill(0.0);
    for (unsigned int d=0;d<dim;d++) refGradU[d]=0.0;

    for (unsigned int k=0; k<dofs_per_cell; ++k){
      if (fe_values.system_to_component_index(k)==DOF){
				for (unsigned int i=0; i<dim; ++i){
	  			if (faceID>-1){
	    			refGradU[i]+=ULocal[k]*fe_values.shape_grad_face(k, q, faceID)[i]; //gradU
	  			}
	  			else{
	    			refGradU[i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //gradU
	  			}
				}
      }
    }
    //Transform gradient to current configuration. gradW=(F^-T)*GradW
    for (unsigned int i=0; i<dim; ++i){
      if (gradientInCurrentConfiguration==false) gradU[q][i]=refGradU[i];
      else{
				gradU[q][i]=0.0;
				for (unsigned int j=0; j<dim; ++j){
	  			gradU[q][i]+=defMap.invF[q][j][i]*refGradU[j];
				}
      }
    }
  }
}

template <class T, int dim>
void model<T, dim>::evaluateVectorFunction(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<2, T>& U,  int faceID){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
 
  //U.fill(0.0); 
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for (unsigned int k=0; k<dim; ++k) U[q][k]=0.0;

    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.system_to_component_index(k) - DOF;
      if (ck>=0 && ck<dim){
				if (faceID>-1){
	 			 U[q][ck]+=ULocal[k]*fe_values.shape_value_face(k, q, faceID); //U
			 	}
			 	else{
				 	U[q][ck]+=ULocal[k]*fe_values.shape_value(k, q); //U
			 	}
     	}
    }
  }
}

template <class T, int dim>
void model<T, dim>::evaluateVectorFunctionGradient(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<3, T>& GradU, int faceID){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  // GradU.fill(0.0);
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for(unsigned int k=0; k<dim; ++k){
      for(unsigned int i=0; i<dim; ++i){
				GradU[q][k][i]=0;
      }
    }
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.system_to_component_index(k) - DOF;
      if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
	  			if (faceID>-1){
	    			GradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad_face(k, q, faceID)[i]; //GradU
	  			}
	 			 	else{
	    		 	GradU[q][ck][i]+=ULocal[k]*fe_values.shape_grad(k, q)[i]; //GradU	    
	    		 // if(!(ULocal[k].val()<100))std::cout<<"ULocal[k]"<<ULocal[k].val()<<std::endl;
	  		 	}
			 	}
      }
    }
  }
}

template <class T, int dim>
void model<T, dim>::evaluateVectorFunctionSecondGradient(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<4, T>& GradGradU, int faceID){
  unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  
  //GradGradU.fill(0.0);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    for(unsigned int k=0; k<dim; ++k){
      for(unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j){
	  			GradGradU[q][k][i][j]=0;
				}
      }
    }
    for (unsigned int k=0; k<dofs_per_cell; ++k){
      unsigned int ck = fe_values.system_to_component_index(k) - DOF;
      if (ck>=0 && ck<dim){
				for (unsigned int i=0; i<dim; ++i){
	  			for (unsigned int j=0; j<dim; ++j){
	   			 	if (faceID>-1){
	      			 GradGradU[q][ck][i][j]+=ULocal[k]*fe_values.shape_grad_grad_face(k, q, faceID)[i][j]; //GradGradU
	    		 	}
	    		 	else{
	      			GradGradU[q][ck][i][j]+=ULocal[k]*fe_values.shape_grad_grad(k, q)[i][j]; //GradGradU
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