#include "../../include/supplementary/dataStruct.h"

template <class T, int dim>
deformationMap<T, dim>::deformationMap(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), detF(n_q_points){}
template deformationMap<Sacado::Fad::DFad<double>, 1>;
template deformationMap<Sacado::Fad::DFad<double>, 2>;
template deformationMap<Sacado::Fad::DFad<double>, 3>;

template <class T, int dim>
deformationMapwithGrad<T, dim>::deformationMapwithGrad(unsigned int n_q_points): F(n_q_points, dim, dim),  invF(n_q_points, dim, dim), gradF(n_q_points, dim, dim, dim), detF(n_q_points){}
template deformationMapwithGrad<Sacado::Fad::DFad<double>, 1>;
template deformationMapwithGrad<Sacado::Fad::DFad<double>, 2>;
template deformationMapwithGrad<Sacado::Fad::DFad<double>, 3>;

