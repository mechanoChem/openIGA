/*
*IGA.h declares class IGA 
*class IGA provides generic commands of FEM at domain level with support of IGAValues and NURBSMesh
*default setting only solve a simple boundary value problem for debugging
*
*for users: Usr should define physics based IGA class derived from this IGA base
*and overwrite virtual functions as necessary  
*/



#ifndef IGA_h
#define IGA_h

//generic headers:deal.ii IGA
#include "headers/headers_dealii.h"
#include "headers/headers_IGA.h"
#include "headers/headers_c.h"
#include <Sacado.hpp>
#include "model/model.h"

using namespace std;

template<int dim>
class IGA
{
public:
	/*
	*IGA constructor and destructor
	*IGA class need NURBS mesh and parameter class as input
	*/
  IGA (NURBSMesh<dim>& _mesh, parametersClass& _params);
  ~IGA();
	
	/*
	*Mark six boundary surfaces
	*/
  void mark_boundaries();
	
	/*
	*Initialize system_rhs, U, Un, dU and system_matrix
	*/
  void setup();
	
	/*
	*cellValue for another way to set fe_values
	*No longer use now, replaced by IGAValues.reinit(*cell)
	*/
  void setupCellValues();
	
	/*
	*Adjust right hand side and Jacobian (system_matrix) to apply dirrchelt BC in solving system
	*/
	void condenseKandRHS();
	
	/*
	*Solve linear system equation
	*direct solver:superLU
	*/
	void solve();
	
	/*
	*top entrance of IGA class 
	*user may overwrite it to initialize their physics based model here
	*/
  virtual void run();
	
	/*
	*apply dirichlet boundary condition
	*user may overwrite it
	*/
  virtual void apply_boundary_conditions();
	
	/*
	*user may overwrite it to apply initial conditions
	*/
  virtual void apply_initial_values();

	/*
	*Fill output vector and output vtk file
	*user can use it directly
	*/
  virtual void output (unsigned int _cycle);
	
	/*
	*Initialization (set system_matrix=0.0; system_rhs=0.0; dU=0.0) before assemble_system_interval
	*Set multiple threads
	*/
  void assemble_system ();
	
	/*
	*for user to assemble system
	*need to overwrite it that uses physics based model model to call getesidual
	*/
  virtual void assemble_system_interval (const typename std::vector<knotSpan<dim> >::iterator &begin, const typename std::vector<knotSpan<dim> >::iterator &end);

	/*
	*pointer to NURBSMesh and model
	*reference of paramtersClass
	*/
  NURBSMesh<dim>* mesh;
  parametersClass& params;
  model<Sacado::Fad::DFad<double>,dim>* modelexample;
	
	/*
	*storage IGAValues for each cell (element)
	*no need to use it now, replaced by IGAValues.reinit(*cell)
	*/
  std::vector<IGAValues<dim>*> cellValues;
	
	/*
	* define sparsity pattern for sparse matrix
	*system_matrix:jacobian materix
	*/
  sparsityPattern sparsity_pattern;
  sparseMatrix system_matrix;
	
	/*
	*system_rhs:right hand side of linear system 
	*U: current solution; Un: solution of previous step iteration; 
	*dU:step change of during iteration used in Newton method
	*/
  denseVector  system_rhs, U, Un, dU;
	
	/*
	*dirichletMap: for applying diirclet boundary condition
	*/
  std::map<unsigned int, double> dirichletMap;

	/*
  *output sequence for output
	*/
  std::vector<solutionClass<dim>* > outputVariables;
	
	/*
	*for loading step by step, numIncrements can be specified NOTE: function run need to overwrite it accordingly
	*currentIncrement: current increment of loading
	*currentIteration:current time step in iteration algorithm 
	*/
	unsigned int numIncrements, currentIncrement, currentIteration;

	/*
	*using deal.iii threads manager
	*/
  dealii::Threads::Mutex     assembler_lock;

};

#endif