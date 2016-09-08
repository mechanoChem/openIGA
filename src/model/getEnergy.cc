#include "../../include/model/model.h"
/*
*for user to access freeEnergy and gradient energy, integration over element
*/
using namespace std;

template <class T, int dim>
double model<T, dim>::getFreeEnergy(){
  return freeEnergy;
}

template <class T, int dim>
double model<T, dim>::getGradientFreeEnergy(){
  return interfaceEnergy;
}

template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;