/*
*IGA<dim>::assemble_system 
*Initialization before assemble_system_interval
*Set multiple threads
*/

#include "../../../include/IGA.h"

using namespace std;

template <int dim>
void IGA<dim>::assemble_system (){
  system_matrix=0.0; system_rhs=0.0; dU=0.0;
  const unsigned int n_threads=dealii::MultithreadInfo::n_threads();
  dealii::Threads::ThreadGroup<> threads;
  typedef typename std::vector<knotSpan<dim> >::iterator knotSpan_iterator;
  std::vector<std::pair<knotSpan_iterator,knotSpan_iterator> > thread_ranges = dealii::Threads::split_range<knotSpan_iterator> (mesh->knotSpanVector.begin(), mesh->knotSpanVector.end(), n_threads);
  printf("start assemble\n");
  for (unsigned int thread=0; thread<n_threads; ++thread){
    threads += dealii::Threads::new_thread (&IGA<dim>::assemble_system_interval, *this, thread_ranges[thread].first, thread_ranges[thread].second);
  }
  threads.join_all ();
  printf("end assemble\n");
  condenseKandRHS();
}

template class IGA<1>;
template class IGA<2>;
template class IGA<3>;




