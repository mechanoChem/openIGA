/**
* Provide sparse data structure
*/


#ifndef SPARSECLASSES_H_
#define SPARSECLASSES_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
#include <map>

#include "../IGAbase/NURBSClass.h"


class sparsityPattern
{	
	private:
		/**
		*initialize with dim=1
		*fill DOFConnections, nzMap, columnIndices and rowIndex
		*/
		void intialize(std::vector<knotSpan<1> >& knotSpanVector);
		
		/**
		*initialize with dim=2
		*fill DOFConnections, nzMap, columnIndices and rowIndex
		*/
		void intialize(std::vector<knotSpan<2> >& knotSpanVector);
		
		/**
		*initialize with dim=3
		*fill DOFConnections, nzMap, columnIndices and rowIndex
		*/
		void intialize(std::vector<knotSpan<3> >& knotSpanVector);
	public:
		/**
		*sparsityPattern class constructor
		*/
		sparsityPattern();
		
		/**
		*initialize with dim=1 for users
		*/
		void init(NURBSMesh<1>* mesh, bool _isMassMatrixPattern=false);
		
		/**
		*initialize with dim=2 for users
		*/
		void init(NURBSMesh<2>* mesh, bool _isMassMatrixPattern=false);
		
		/**
		*initialize with dim=3 for users
		*/
		void init(NURBSMesh<3>* mesh, bool _isMassMatrixPattern=false);
		
		/**
		*map of non-zeros/zero entries
		*/
		std::vector<std::map<unsigned int, unsigned int> > nzMap;
		
		/**
		*column Indices
		*/
		std::vector<int> columnIndices;
		
		/**
		*row indices
		*/
		std::vector<int> rowIndex;
		/**
		*nnz:nonzero value size
		*/
		unsigned int nnz; bool isMassMatrixPattern;
		std::map<unsigned int, std::set<unsigned int> > DOFConnections;
};

/**
*denseMatrix class
*/
class denseMatrix
{	
	public:
		denseMatrix(unsigned int _size1, unsigned int _size2);
		
		/**
		*values store all value of the denseMatrix
		*/
		std::vector<std::vector<double> > values;
		
		/**
		*size of this denseMatrix
		*/
		unsigned int size();
		
		/**
		*operator  overloading
		*(i,j) access the reference of value at row i and column j
		*/
		double& operator()(unsigned int _a, unsigned int _b);
		
		/**
		*set denseMatrix to be zero that is zero number with value "0"
		*/
		void operator= (double _value);
		
		/**
		*print out all value of the denseMatrix in matrix formula
		*/
		void print();
};

/**
*denseVector class
*/
class denseVector
{	
	public:
		/**
		*constructors
		*/
		denseVector();
		
		/**
		*constructors with _size
		*/
		denseVector(unsigned int _size);
		
		/**
		*values store all value of the denseVector
		*/
		std::vector<double> values;
		
		/**
		*size of the vector
		*/
		unsigned int size();
		
		/**
		*L2_norm: sum of (a_i)^2
		*/
		double l2_norm();
		
		/**
		*generate denseVector corresponding to sparse Matrix
		*/
		void reinit(sparsityPattern& sparsitypattern);
		
		/**
		*generate dsenseVector with size "_size"
		*/
		void reinit(unsigned int _size);
		
		/**
		*return reference of value at _a
		*/
		double& operator() (unsigned int _a);
		
		/**
		*set denseVector to be zero that is zero number with value "0"
		*/
		void operator= (double _value);
		
		/**
		*set two senseVector to be same
		*/
		void operator= (denseVector& _value);
		
		/**
		*add denseVector (input) to this denseVector
		*/
		void operator+= (denseVector& _value);
		
		/**
		*print out all value of the denseVector
		*/
		void print();
};

/**
*sparseMatrix class
*/
class sparseMatrix
{
	private:
		sparsityPattern* sparsitypattern;
	public:
		std::vector<double> nonZeroValues;
		/**
		* class constructors
		*/
		sparseMatrix();
		void reinit(sparsityPattern& _sparsitypattern);
		
		/**
		*operator  overloading
		*(i,j) access the reference of value at row i and column j
		*/
		double& operator() (unsigned int _a, unsigned int _b);
		
		/**
		*return value at row i and column j
		*/
		double val (unsigned int _a, unsigned int _b);
		
		/**
		*set sparseMatrix to be zero that is zero number with value "0"
		*/
		void operator= (double _value);
		
		/**
		*print out all value of the sparseMatrix in matrix formula
		*/
		void print();
};


#endif