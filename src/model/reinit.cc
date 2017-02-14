#include "../../include/model/model.h"
/*
*re-initialize parameters for class model
*/
template <class T, int dim>
void model<T, dim>::reinit(parametersClass<dim>& _params)
{
	params=&_params;
  lambda= params->getDouble("lambda");
  mu=  params->getDouble("mu");
  muSG= params->getDouble("muSG");
  gamma= params->getDouble("Gamma");
  C= params->getDouble("C");
	NumKnotInterval=params->getInt("knots");
  he=1.0/NumKnotInterval;
	DOF=params->getInt("DOF");
	finiteStrain=params->getBool("finiteStrain");
	
  bcType= params->getString("bcType").c_str();
}
template class model<Sacado::Fad::DFad<double>,1>;
template class model<Sacado::Fad::DFad<double>,2>;
template class model<Sacado::Fad::DFad<double>,3>;