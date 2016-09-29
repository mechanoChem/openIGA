/**
*Model.h declares base model class.
*class model provides generic quadrature level commands. 
*		1. GetResidual for mechanics and B.C.
*		2. Evaluate deformation and gradient
*functions could be overloeaded
*
*for user:user should define their physics based model derived from this base model. 
*/


#ifndef MODEL_h
#define MODEL_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
//generic headers:deal.ii, IGA
#include "deal.II/base/table.h"
#include <Sacado.hpp>
#include "../supplementary/dataStruct.h"
#include "../supplementary/parameters.h"
#include "../IGAbase/NURBSClass.h"

using namespace std;

template <class T, int dim>
class model
{
public:
	
 	/**
 	*model constructor
 	*/
  model ();
	
 	/**
 	*model destructor
 	*/
  ~model();
	
	/**
	*re-initialize parameter class
	*/
	void reinit(parametersClass& _params);
	
	/**
	*main entrance to call all individual "residual" functions
	*user may overwrite it
	*/
  virtual void getResidual(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R,unsigned int currentIteration);

	/**
	*residual for regular stress and high order stress
	*/
	void residualForMechanics(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R);

	/**
	*residual for neumman boundary condition
	*user should overwrite it
	*/
	virtual void residualForNeummanBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R);
	
	/**
	*residual for high order boundary condition
	*User may use it directly though it is virtual
	*/
	virtual void residualForHighOrderBC(knotSpan<dim>& cell, IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<1, T >& R);
	
	/**
	*Evaluate deformation gradient at given quadrature points
	*Calculate Piola- kirchho stress and higher order stress Beta
	*/
	virtual void evaluateStress(IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, dealii::Table<3, T>& P, dealii::Table<4, T>& Beta,int faceID=-1);
	
	/**
	*Evaluate gradient of deformation tensor
	*Evaluate determinate and inverse of deformation tensor
	*/
	void getDeformationMapWithGradient(IGAValues<dim>& fe_values, dealii::Table<1, T >& ULocal, deformationMapwithGrad<T, dim>& defMap, int faceID=-1);
	
	/**
	*Evaluate scalar functions using quadrature points (first three dof e.g. displacements) interpolation
	*/
	void evaluateScalarFunction(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<1, T>& U,  int faceID =-1);
	
	/**
	*Evaluate scalar functions gradient using quadrature points (first three dof e.g. displacements) interpolation 
	*gradientInCurrentConfiguration is to be actived when evaluate gradient at current configuration
	*/
	void evaluateScalarFunctionGradient(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<2, T>& gradU, deformationMap<T, dim>& defMap, bool gradientInCurrentConfiguration, int faceID =-1);
	
	/**
	*Evaluate vector functions using quadrature points (first three dof e.g. displacements) interpolation
	*/
	void evaluateVectorFunction(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<2, T>& U,  int faceID =-1);
	
	/**
	*Evaluate vector functions gradient using quadrature points (first three dof e.g. displacements) interpolation
	*/ 
	void evaluateVectorFunctionGradient(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<3, T>& GradU, int faceID =-1);
	
	/**
	*Evaluate vector functions second gradient using quadrature points (first three dof e.g. displacements) interpolation
	*/ 
	void evaluateVectorFunctionSecondGradient(IGAValues<dim>& fe_values, dealii::Table<1, T>& ULocal, dealii::Table<4, T>& GradGradU, int faceID =-1);
	
	/**
	*return freeEnergy integration over element
	*/
	double getFreeEnergy();
	
	/**
	*return gradient energy integration over element
	*/
	double getGradientFreeEnergy();
	
	/**
	*pointer to material paramters
	*/
	parametersClass* params;
	
	/**
	*current dof of residual
	*DOF=0 (default) for evaluating residual of mechanical term  
	*/
	unsigned int DOF;
	
	/**
	*finite strain=true, using finite strain elasticity
	*/
	bool finiteStrain;
	
	/**
	*current iteration
	*/
	unsigned int iteration;
	
	/**
	*freeEnergy integration over element
	*/
	double freeEnergy;
	
	/**
	*gradient energy integration over element
	*/
	double interfaceEnergy;
	
	/**
	*material parameters

	*deatials explanations of these parameters can be found at "Rudraraju, Van der Ven, Garikipati, Comp. Meth. App. Mech. Engrg., 278 705, 2014,
	*title: Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains" 
	*/
  double lambda;
	
	/**
	*material parameters
	*elasticity constant
	*deatials explanations of these parameters can be found at "Rudraraju, Van der Ven, Garikipati, Comp. Meth. App. Mech. Engrg., 278 705, 2014,
	*title: Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains" 
	*/
  double mu;
	
	/**
	*material parameters
	*muSG:mu for strain gradient term(beta) usually is same with mu
	*deatials explanations of these parameters can be found at "Rudraraju, Van der Ven, Garikipati, Comp. Meth. App. Mech. Engrg., 278 705, 2014,
	*title: Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains" 
	*/
  double muSG;
	
	/**
	*material parameters
	*l is gradient length scale
	*deatials explanations of these parameters can be found at "Rudraraju, Van der Ven, Garikipati, Comp. Meth. App. Mech. Engrg., 278 705, 2014,
	*title: Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains" 
	*/
	double l;
	
	/**
	*material parameters
	*gamma is parameter for weakly enforces the higher-order Dirichlet boundary condition using a penalty-based approach
	*deatials explanations of these parameters can be found at "Rudraraju, Van der Ven, Garikipati, Comp. Meth. App. Mech. Engrg., 278 705, 2014,
	*title: Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains" 
	*/
  double gamma;
	
	/**
	*material parameters
	*C is parameter for weakly enforces the higher-order Dirichlet boundary condition using a penalty-based approach
	*deatials explanations of these parameters can be found at "Rudraraju, Van der Ven, Garikipati, Comp. Meth. App. Mech. Engrg., 278 705, 2014,
	*title: Three dimensional iso-geometric solutions to general boundary value problems of Toupin's theory of gradient elasticity at finite strains" 
	*/
  double C;
	
	/**
	*NumKnotInterval:number of knot interval used in weakly applying higher order boundary condition
	*/
	double NumKnotInterval;
	
	/**
	*he=1.0/NumKnotInterval;
	*/
  double he;
	
	/**
	*may only used for example (redundant)
	*/ 
  double load;
	
	/**
	*type of boundary condition
	*/ 
  const char* bcType;

};

#endif
