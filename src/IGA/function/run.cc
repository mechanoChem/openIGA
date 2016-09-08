/*
*IGA<dim>::run
*top entrance of IGA
*for user to initialize their physics based model
*/

#include "../../../include/IGA.h"
#include <Sacado.hpp>

using namespace std;

template<int dim>
void IGA<dim>::run (){
	numIncrements=1; currentIncrement=0;
  model<Sacado::Fad::DFad<double>,dim> _Modelexample;
	modelexample=&_Modelexample;
	modelexample->reinit(params);
	
  setup();
  mark_boundaries();
	apply_initial_values();
	apply_boundary_conditions();
	solve(); 
	output(0);
}

template class IGA<1>;
template class IGA<2>;
template class IGA<3>;