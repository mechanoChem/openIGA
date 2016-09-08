#include "../../include/model/model.h"
/*
*residual for high order boundary condition
*/
template <class T, int dim>

void model<T, dim>::residualForHighOrderBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R)
{
  //Weak dirichlet condition grad(u).n=0
	unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  for (unsigned int faceID=0; faceID<2*dim; faceID++){
		if ((cell.boundaryFlags[faceID]==(dim-1)*2+1) or (cell.boundaryFlags[faceID]==(dim-1)*2+2)){
		//compute face normal (logic for computing n like this only works for cube geometries)
		std::vector<double> n(dim, 0);
		n[(cell.boundaryFlags[faceID]-1)/2]=std::pow(-1.0, (int)cell.boundaryFlags[faceID]%2);
	
		//Temporary arrays
		dealii::Table<3,Sacado::Fad::DFad<double> > PFace (fe_values.n_face_quadrature_points, dim, dim);
		dealii::Table<4,Sacado::Fad::DFad<double> > BetaFace (fe_values.n_face_quadrature_points, dim, dim, dim);
		evaluateStress(fe_values, ULocal, PFace, BetaFace,faceID);
 	
		//evaluate gradients on the faces
		dealii::Table<3,Sacado::Fad::DFad<double> > uij(fe_values.n_face_quadrature_points, dim, dim); //uij.fill(0.0);
		dealii::Table<4,Sacado::Fad::DFad<double> > uijk(fe_values.n_face_quadrature_points, dim, dim, dim); //uijk.fill(0.0);
		for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
      for (unsigned int i=0; i<dim; ++i){
        for (unsigned int j=0; j<dim; ++j){
          uij[q][i][j]=0.0;
          for (unsigned int k=0; k<dim; ++k){
            uijk[q][i][j][k]=0.0;
           }
         }
       }
     }
		 //Loop over quadrature points
		 for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
			 for (unsigned int d=0; d<dofs_per_cell; ++d){
				 unsigned int i = fe_values.system_to_component_index(d) - DOF;
	    	 for (unsigned int j=0; j<dim; ++j){
	      	 uij[q][i][j]+=ULocal[d]*fe_values.shape_grad_face(d, q, faceID)[j];
	      	 for (unsigned int k=0; k<dim; ++k){
						 uijk[q][i][j][k]+=ULocal[d]*fe_values.shape_grad_grad_face(d, q, faceID)[j][k];
					 }
	    	  }
	  	  }	
			}

			//evaluate tensor multiplications 
			dealii::Table<2,Sacado::Fad::DFad<double> > uijn(fe_values.n_face_quadrature_points, dim); //uijn.fill(0.0);
			dealii::Table<2,Sacado::Fad::DFad<double> > Betaijknn(fe_values.n_face_quadrature_points, dim);// Betaijknn.fill(0.0);
			dealii::Table<2,double> wijn(fe_values.n_face_quadrature_points, dofs_per_cell); //wijn.fill(0.0);
			dealii::Table<2,double> wijknn(fe_values.n_face_quadrature_points, dofs_per_cell); //wijknn.fill(0.0);
      for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
        for (unsigned int i=0; i<dim; ++i){
          uijn[q][i]=0.0;
          Betaijknn[q][i]=0.0;
         }
         for(unsigned int j=0;j<dofs_per_cell;j++){
           wijn[q][j]=0.0;
           wijknn[q][j]=0.0;
         }
       }
			//evaluate uijn, Betaijknn
			for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  		double tempdirchletBC=0.0;
	  		for (unsigned int i=0; i<dim; ++i){
	    		for (unsigned int j=0; j<dim; ++j){
	      		uijn[q][i]+=uij[q][i][j]*n[j];
	      		for (unsigned int k=0; k<dim; ++k){
							Betaijknn[q][i]+=BetaFace[q][i][j][k]*n[j]*n[k];
	      		}
	    		}
	    		//tempdirchletBC+=std::pow(uijn[q][i].val(),2.0);
	  		}
	  		//dirchletBC=std::max(dirchletBC, std::pow(tempdirchletBC,0.5));
	  		//gradun(fe_values.cell->id,q,i)= uijn[q][i].val();
			 }
			 //evaluate wijn, wijknn
			 for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){
	  		 for (unsigned int d=0; d<dofs_per_cell; ++d){
	    		 for (unsigned int j=0; j<dim; ++j){
	      		 wijn[q][d]+= fe_values.shape_grad_face(d, q, faceID)[j]*n[j];
	      			for (unsigned int k=0; k<dim; ++k){
								wijknn[q][d]+= fe_values.shape_grad_grad_face(d, q, faceID)[j][k]*n[j]*n[k];
	      			}
	    			}
	  			}	  
				}
	
				//Add the weak dirichlet terms
				for (unsigned int d=0; d<dofs_per_cell; ++d) {
	  			unsigned int i = fe_values.system_to_component_index(d) - DOF;
	  			if (i==dim-1){ //enforcing weak dirichlet only along dim direction
	    		  for (unsigned int q=0; q<fe_values.n_face_quadrature_points; ++q){	
	      			//-mu*l*l*w_(ij)n_j*Beta_(ijk)n_jn_k
	      			R[d] += -wijn[q][d]*Betaijknn[q][i]*fe_values.JxW_face(q, faceID);
	      			//-mu*l*l*gamma*wijknn*uijn //WRONG form curretnly. Not used now as there is no point trying to make an unsymmetric euqation adjoint consistent. So gamma is always set to zero.
	      			gamma=0.0;
	      			R[d] += -gamma*wijknn[q][d]*uijn[q][i]*fe_values.JxW_face(q, faceID);
	      			//mu*l*l*(C/he)*wijn*uijn
	      			R[d] += (C/he)*wijn[q][d]*uijn[q][i]*fe_values.JxW_face(q, faceID);
	    			}
	  			}
				} 
    }
  } 
}

template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;