#include "../../include/model/model.h"
/*
*Calculate deformation gradient at given quadrature points
*Evaluate Piola- kirchho stress and higher order stress Beta
*/
template <class T, int dim>
void model<T, dim>::evaluateStress(IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<3, T>& P, dealii::Table<4, T>& Beta, int faceID){
  //unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;

  deformationMapwithGrad<Sacado::Fad::DFad<double>, dim> defMap(n_q_points); 
  getDeformationMapWithGradient(fe_values,ULocal, defMap, faceID);
  
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){    
    //Compute F
    T F[dim][dim], dF[dim][dim][dim];
    for (unsigned int i=0; i<dim; i++) {
      for (unsigned int j=0; j<dim; j++) {
				F[i][j]=defMap.F[q][i][j];	
				for (unsigned int k=0; k<dim; k++) {
	  			dF[i][j][k]=defMap.gradF[q][i][j][k];
				}
      }
    }
    //Compute strain metric, E  (E=0.5*(F^T*F-I))
    T E[dim][dim], dE[dim][dim][dim], trE=0.0, trE2=0.0;
    for (unsigned int I=0; I<dim; I++){
      for (unsigned int J=0; J<dim; J++){
				E[I][J] = -0.5*(I==J);
				for (unsigned int k=0; k<dim; k++){
	  			E[I][J] += 0.5*F[k][I]*F[k][J];
				}
				for (unsigned int K=0; K<dim; K++){
	  			dE[I][J][K]=0.0;
	  			for (unsigned int k=0; k<dim; k++){
	    			dE[I][J][K] += 0.5*(F[k][I]*dF[k][J][K]+F[k][J]*dF[k][I][K]);
	  			}
				}
				trE2 += E[I][J]* E[I][J];	
      }
      trE +=E[I][I];   
    }

    //infinitesimal strain implementation
    if (!finiteStrain){
      trE=0.0; trE2=0.0;
      for (unsigned int i=0; i<dim; i++) {
				for (unsigned int j=0; j<dim; j++) {
	  			F[i][j]=(i==j);	
	  			E[i][j]=(defMap.F[q][i][j] - (i==j)+defMap.F[q][j][i] - (i==j))/2;
	  			for (unsigned int k=0; k<dim; k++) {
	    			dE[i][j][k]=(defMap.gradF[q][i][j][k]+defMap.gradF[q][j][i][k])/2;
	  			}
	  			trE2 += E[i][j]*E[i][j];	
				}
				trE +=E[i][i];   
      }
    }

    //compute P and Beta and energy
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int J=0; J<dim; ++J){
				P[q][i][J]=lambda*trE*F[i][J];
				for (unsigned int K=0; K<dim; ++K){
	  			P[q][i][J]+= 2*mu*F[i][K]*E[K][J];
	  			Beta[q][i][J][K]=0.0;
	  			for (unsigned int A=0; A<dim; ++A){
	    			P[q][i][J]+= mu*l*l*dF[i][K][A]*dE[J][K][A];
	    			Beta[q][i][J][K]+=muSG*l*l*F[i][A]*dE[A][J][K];
	  			}
				}
      }
    }
		
		//calculate Energy
    if (faceID==-1){
      double tempFreeEnergy=0.0, tempInterfaceEnergy=0.0;
      tempFreeEnergy= 0.5*lambda*std::pow(trE.val(),2)+ mu*trE2.val();
      for (unsigned int I=0; I<dim; ++I){
				for (unsigned int J=0; J<dim; ++J){
	  			for (unsigned int K=0; K<dim; ++K){
	    			tempInterfaceEnergy+=0.5*muSG*l*l*std::pow(dE[I][J][K].val(),2);
	  			}	
				}
      }
      freeEnergy+=tempFreeEnergy*fe_values.JxW(q);
      interfaceEnergy+=tempInterfaceEnergy*fe_values.JxW(q);
    }
  }
}

template class model<Sacado::Fad::DFad<double>, 1>;
template class model<Sacado::Fad::DFad<double>, 2>;
template class model<Sacado::Fad::DFad<double>, 3>;