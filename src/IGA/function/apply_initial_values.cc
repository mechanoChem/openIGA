/*
*IGA<dim>::apply_initial_values
*for user to apply initial conditions
*Un is solution of previous step iteration; U is current solution
*/
 
 #include "../../../include/IGA.h"
using namespace std;

//Apply initial conditions
template<int dim>
void IGA<dim>::apply_initial_values(){
  Un=0.0; U=0.0;
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;