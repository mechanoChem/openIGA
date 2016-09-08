#include "../../include/model/model.h"
/*
*Evaluate gradient of deformation tensor
*Evaluate determinate and inverse of deformation tensor
*/
template <class T, int dim>
void model<T, dim>::getDeformationMapWithGradient(IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, deformationMapwithGrad<T, dim>& defMap, int faceID)
{
  //unsigned int dofs_per_cell= fe_values.dofs_per_cell;
  unsigned int n_q_points= (faceID>-1) ?  fe_values.n_face_quadrature_points : fe_values.n_quadrature_points;
  //evaluate dx/dX
  dealii::Table<3, T> GradU(n_q_points, dim, dim);
  dealii::Table<4, T> GradGradU(n_q_points, dim, dim, dim);
  evaluateVectorFunctionGradient(fe_values, ULocal, GradU, faceID);
  evaluateVectorFunctionSecondGradient(fe_values, ULocal, GradGradU, faceID);
  //Loop over quadrature points
  for (unsigned int q=0; q<n_q_points; ++q){
    dealii::Table<2, T > Fq(dim, dim), invFq(dim, dim); T detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.F[q][i][j] = Fq[i][j] = (i==j) + GradU[q][i][j]; //F (as double value)
				for (unsigned int k=0; k<dim; ++k){
	  			defMap.gradF[q][i][j][k] = GradGradU[q][i][j][k]; //gradF 
				}
      }
    }
    getInverse<T, dim>(Fq, invFq, detFq); //get inverse(F)
    defMap.detF[q] = detFq;
    for (unsigned int i=0; i<dim; ++i){
      for (unsigned int j=0; j<dim; ++j){
				defMap.invF[q][i][j] = invFq[i][j];
      }
    }
    //detF
    if (defMap.detF[q].val()<=1.0e-15 && iteration==0){
      std::vector<double> coords(fe_values.quadraturePoints_coords(q));
      printf("**************Non positive jacobian detected**************. Value: %12.4e\n", defMap.detF[q].val());
      for (unsigned int i=0; i<dim; ++i){
				for (unsigned int j=0; j<dim; ++j) printf("%12.6e  ", defMap.F[q][i][j].val());
					printf("\n"); exit(-1);
      	}
      throw "Non positive jacobian detected";
    }
  }  
}
	
template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;