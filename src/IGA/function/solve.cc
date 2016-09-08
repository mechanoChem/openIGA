/*
*IGA<dim>::solve
*Solve linear system equation
*direct solver
*/

#include "../../../include/IGA.h"

using namespace std;

template<int dim>
void IGA<dim>::solve (){
  double res=params.getDouble("res");
	double tol=params.getDouble("tol");
	double abs_tol=params.getDouble("abs_tol");
	double initial_norm=params.getDouble("initial_norm");
	double current_norm=params.getDouble("current_norm");
	double max_iteration=params.getDouble("max_iteration");
  currentIteration=0;
  double oldNorm=1;
  while (true){
    if (currentIteration>=max_iteration){printf ("Maximum number(default 1) of iterations reached without convergence. \n"); break; exit (1);}
    double currentNorm=dU.l2_norm();
    if (current_norm>1/std::pow(tol,2)){printf ("\nNorm is too high. \n\n"); exit (1);}
    assemble_system();
    current_norm=system_rhs.l2_norm(); 
    initial_norm=std::max(initial_norm, current_norm);
    res=current_norm/initial_norm;
    printf ("Inc:%3u, Iter:%2u. Residual norm: %10.2e. Relative norm: %10.2e, UNorm: %12.8e \n", currentIncrement, currentIteration, current_norm, res, currentNorm); 
    double finalNorm=res;
    if (res<tol || current_norm< abs_tol){printf ("Residual converged in %u iterations.\n", currentIteration); break;}
		std::cout<<"before solve"<<std::endl;
     int status=luSolver(&system_matrix.nonZeroValues.at(0), &sparsity_pattern.columnIndices.at(0), &sparsity_pattern.rowIndex.at(0), &system_rhs.values.at(0), (int) system_rhs.size(), (int) system_rhs.size(), (int) sparsity_pattern.nnz, &dU.values.at(0), 0);
    U+=dU;
    ++currentIteration;
  }
  Un=U;	
}
template class IGA<1>;
template class IGA<2>;
template class IGA<3>;