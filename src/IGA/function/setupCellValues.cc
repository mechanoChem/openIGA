/*
*IGA<dim>::setupCellValues
*Set cellValue for another way to set fe_values
*NO longer use now
*/

#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::setupCellValues(){
  printf("setting up cell values\n");
  for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
    IGAValues<dim>* tempCellValues=new IGAValues<dim>(mesh, dim, 2);
    tempCellValues->reinit(*cell);
    cellValues.push_back(tempCellValues);
  }
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;