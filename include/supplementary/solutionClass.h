/*
* Provide storage and projection for FEM results with NURBS mesh
*data type include scalar, vector and tensor
*/

#ifndef SOLUTIONCLASSES_H_
#define SOLUTIONCLASSES_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>

#include "../headers/headers_solver.h"

#include "sparsityDataClass.h"
#include "../IGAbase/NURBSClass.h"

using namespace std;

enum dataType {SCALAR, VECTOR, TENSOR};
enum dataLocation {NODAL, QUADRATURE};

template <int dim>
class solutionClass{
	private:
		NURBSMesh<dim>* mesh;
		
	public:
	 	/*
	 	*class constructor and destructor
		*resize values vector according to data type
	 	*/
		solutionClass(NURBSMesh<dim>& _mesh, dataLocation location, dataType data, std::string  _variableName);
		~solutionClass();
		
		/*
		*create mass jacobian
		*/
		void createMassMatrix();
		
		/*
		*project result onto nodes from quadrature point using quarature points interpolation
		*/
		void projectQuadratureValues();
		
		/*
		*vector to store initial data value (at nodes or quadrature points)
		*/
		std::vector<double> values;
		
		/*
		*vector to store sata value after projected onto nodes
		*/
		std::vector<double> projectedValues;
		
		/*
		*dataLocation: at NODAL or QUADRATURE
		*datatype:SCALAR, VECTOR or TENSOR
		*/
		dataLocation datalocation;
		dataType datatype;
		
		/*
		*operator  overloading
		*return values 
		*/
		double& operator() (unsigned int _a) {
			if ((datalocation==NODAL) && (datatype==SCALAR)) return values.at(_a);
			else {printf("solutionClass: incompatible arguments\n 1"); exit(-1);}
		}
		double& operator() (unsigned int _a, unsigned int _b) {
			if (datalocation==NODAL){
				if (datatype==VECTOR) return values.at(_a*dim+_b);
				else if (datatype==TENSOR) return values.at(_a*dim*dim+_b);
				else {printf("solutionClass: incompatible datatype\n"); exit(-1);}
			}
			else if ((datalocation==QUADRATURE) && (datatype==SCALAR)) return values.at(_a*numQuadPoints +_b);
			else {printf("solutionClass: incompatible arguments\n 2"); exit(-1);}
		}
		double& operator() (unsigned int _a, unsigned int _b, unsigned int _c) {
			if (datalocation==QUADRATURE){
				if (datatype==VECTOR) return values.at(_a*numQuadPoints*dim+_b*dim+_c);
				else if (datatype==TENSOR) return values.at(_a*numQuadPoints*dim*dim+_b*dim*dim+_c);
				else {printf("solutionClass: incompatible datatype\n 3"); exit(-1);}
			}
			else {printf("solutionClass: incompatible arguments\n 4"); exit(-1);}
		}
		
		/*
		*numQuadPoints=mesh->quadPtStencilSize 
		*numVariablesPerPoint=dim for vector;=dim*dim for tensor
		*/
		unsigned int numQuadPoints;
		unsigned int numVariablesPerPoint;
		/*
		*solution name to be shown in vtk file
		*/
		std::string variableName;
		
		/*
		*massMatrix and corresponding sparsityPattern 
		*only for projection use
		*/
		static sparseMatrix* massMatrix;
		static sparsityPattern* mass_sparsity_pattern;
};

#endif 