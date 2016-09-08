/*
*IGA<dim>::mark_boundaries
*Mark six boundary surfaces
*/

#include "../../../include/IGA.h"

using namespace std;

//Mark boundary flags
template<int dim>
void IGA<dim>::mark_boundaries(){
  for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
    for (unsigned int i=0;i<dim; i++){
      if (cell->endKnots[i][0]==0) cell->boundaryFlags[i*2+0]=i*2+1;
      if (cell->endKnots[i][1]==1) cell->boundaryFlags[i*2+1]=i*2+2;
    }
  }
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;