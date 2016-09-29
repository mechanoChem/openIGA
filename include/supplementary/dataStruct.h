/**
* Provide physics based data structure
*/

#ifndef FUNCTIONEVALUATIONS_H_
#define FUNCTIONEVALUATIONS_H_
#include "deal.II/base/table.h"

template <class T, int dim>
  struct deformationMap{
  deformationMap(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), detF(n_q_points){}
    dealii::Table<3, T> F, invF;
    dealii::Table<1, T> detF;
 };

template <class T, int dim>
  struct deformationMapwithGrad{
  deformationMapwithGrad(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), gradF(n_q_points, dim, dim, dim), detF(n_q_points){}
    dealii::Table<3, T> F, invF;
    dealii::Table<4, T> gradF;
    dealii::Table<1, T> detF;
  };
	

#endif	