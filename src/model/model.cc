#include "../../include/model/model.h"

template <class T, int dim>
model<T, dim>::model ():freeEnergy(0),interfaceEnergy(0)
{
	std::cout<<"MODEL GENERATED"<<std::endl;
}

template <class T, int dim>
model<T, dim>::~model (){}


template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;
