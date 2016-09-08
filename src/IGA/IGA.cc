#include "../../include/IGA.h"

using namespace std;

template<int dim>
IGA<dim>::IGA (NURBSMesh<dim>& _mesh,parametersClass& _params): mesh(&_mesh),params(_params)
{

}

template<int dim>
IGA<dim>::~IGA (){}

template class IGA<1>;
template class IGA<2>;
template class IGA<3>;