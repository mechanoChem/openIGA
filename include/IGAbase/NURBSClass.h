/**
*NURBSCLASS.h declares all IGAbase classes which are depended on each other.
*class nativeNURBSStructure reads .h5 and collects geometry information
*
*class NURBSMesh reads geometry and generates NURBS mesh (elements or knotspan). 
*  On element(knotspan) level we have
*			1. lineflags to support mark element edge
*			2. planflags to support mark plane inside element
*			3. defectflags to support mark element contain certain defects (for our defects problem)
*
*class IGAValues
*  1. reads NURBS mesh
*  2. calculates NURBS basis functions
*  3. provides generic FEM commands.
*
*class basisFunction, controlPoint and knotspan provides support for NURBSMesh and IGAValurs
*For user: Nothing should be touched here
*/


#ifndef NURBSCLASSES_H_
#define NURBSCLASSES_H_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include <vector>
#include <set>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/tensor.h>
#include <H5Cpp.h>
#include "../supplementary/supplementaryFunctions.h"

//forward declaration
class sparsityPattern;

template <int dim>
class IGAValues;

template <int dim>
class NURBSMesh;

template <int dim>
class nativeNURBSStructure
{
	public:
		/**
		*nativeNURBSStructure constructor
		*nativeNURBSStructure class need .h5 file as input
		*/
		nativeNURBSStructure(const char* _hdf5FileName);
		
		/**
		*nativeNURBSStructure destructor
		*/
		~nativeNURBSStructure();
		
		/**
		*read .h5 mesh and fill mesh information
		*/
		int readHDF5();
		
		/**
		*hdf5FileName=_hdf5FileName
		*/
		const char* hdf5FileName;
		
		/**
		*polynomial order along each direction
		*/
		std::vector<unsigned int> order;
		
		/**
		*number of control points along each direction
		*/
		std::vector<unsigned int> numCtrlPts;
		
		/**
		knot values along each direction
		*/
		std::vector<std::vector<double> > knots;
		
		/**
		*control point coefficients
		*/
		std::vector<std::vector<double> > ctrlPtCoefs;
		
		/**
		*central point of knot along each direction
		*/
		std::vector<std::vector<unsigned int> > ctrlPtKnots;
};

/**
*NURBS basis function formula
*/
template <int dim>
class basisFunction
{
	private:
		/**
		*recursive evaluation of value of NURBS using Cox-DeBoor recursion formula
		*/
		double evaluateN(unsigned int i, unsigned int localKnotIndex, unsigned int _p, double xi, unsigned int k=0);
	public:
		
		/**
		*basisFunction constructor and destructor
		*/
		basisFunction();
		
		/**
		*weight of control point
		*weight=mesh.ctrlPtCoefs[dim].at(id);
		*/
		double weight;
		
		/**
		*p_i: order  along each i^th direction
		*p[a]=mesh.order[a]
		*/
		std::vector<unsigned int> p; 
		
		/**
		*p_i+2 localKnots vector along each i^th direction
		*store coordinates of control points inside local knots
		*/
		std::vector<std::vector<double> > localKnots; 

		
		/**
		*detect control point which belongs to the knotspan e.g. value is nonzero
		*/
		bool pointInBasisSupport(std::vector<double>& coords);
		
		/**
		*detect control point inside knot span (danger in use)
		*/
		bool pointCrossSurfaceSupport(double x);
		
		/**
		*evaluate value of basisfunction:N*weight
		*/
    double value(std::vector<double>& coords, std::vector<unsigned int>& dervs);
};

/**
*controlPoint class
*/
template <int dim>
class controlPoint
{
	public:
		/**
		*controlPoint constructor
		*_id is number of control point
		*/
		controlPoint(unsigned int _id, nativeNURBSStructure<dim>& mesh, std::vector<unsigned int>& _dofs);
		
		
		unsigned int id;
		/**
		*dofPerControlPoint defined by users
		*/
		std::vector<unsigned int> dofs;
		
		/**
		*coordinate of point point
		*/
		std::vector<double> coords;
		/**
		*each control point has a basis class
		*/
		basisFunction<dim> basis; 
};

/**
*knotspan class
*/
template <int dim>
class knotSpan
{
	public:
		/**
		*knotSpan constructor
		*/
		knotSpan(unsigned int _id,  nativeNURBSStructure<dim>& mesh, NURBSMesh<dim>& nurbsmesh, std::vector<unsigned int>& tempIndices, IGAValues<dim>& fe_values); 
		
		/**
		*id:number of control point
		*/
		unsigned int id;
		
		/**
		*boundaryFlags: denote cells which contain boundary surface
		*/
		std::vector<unsigned int> boundaryFlags;
		
		/**
		*lineFlags: mark line of cell
		*/
		std::vector<unsigned int> lineFlags;
		
		/**
		*planeFlags mark plane inside of cell
		*/
		std::vector<unsigned int> planeFlags;
		
		/**
		*defectFlags: mark cells that contain defects 
		*/
		std::vector<unsigned int> defectFlags;
		
		/**
		*for manually specify one special quadrature point coordinate along one direction of plane
		*/
		double Sp_planeQuad_Point;
		
		/**
		*for manually specify special quadrature point along two direction on the line
		*/
    std::vector<double> Sp_lineQuad_Point;
		
		/**
		*coordinates of end points of knots
		*/
		std::vector<std::vector<double> > endKnots;
		
		/**
		*coordinate of edge on surfaces of knots
		*/
		std::vector<std::vector<double> > edgeCoords;
		
		/**
		*coordinate of center on surfaces of knots
		*/
		std::vector<double> center;
		
		
    std::vector<double> controlPoint_cross_surface_cof;
		
		/**
		*control points that in basis support (non-zero) of the knot
		*/
		std::vector<controlPoint<dim>* > controlPoints;
		
		/**
		*generate local_dof_indices to global_dof_indices
		*/
		std::vector<unsigned int> local_dof_indices;
		
		/**
		*generate local_mass_dof_indices (local indices of control point) to global_dof_indices 
		*i.e. same with local_dof_indices if each control point only have one dof
		*/
		std::vector<unsigned int> local_mass_dof_indices;
		
		/**
		*determine wether points inside the knot span
		*/
		bool pointInKnotSpan(std::vector<double>& coords);
};

/**
*NURBSMesh class
*/
template <int dim>
class NURBSMesh
{
	public:
		/**
		*NURBSMesh constructor
		*/
		NURBSMesh(nativeNURBSStructure<dim>& mesh, unsigned int _dofPerControlPoint, unsigned int _quadPtStencilSize);
		
		/**
		*NURBSMesh destructor
		*/
		~NURBSMesh();
		
		/**
		*dofPerControlPoint and quadPtStencilSize are defiend by users
		*/
		unsigned int dofPerControlPoint, quadPtStencilSize;
		
		/**
		*NURBSMesh contains all information of each control point
		*/
		std::vector<controlPoint<dim> > controlPointVector;
		
		/**
		*NURBSMesh contains all information of each knot span 
		*/
		std::vector<knotSpan<dim> > knotSpanVector;
};

/**
*IGA value class
*/
template <int dim>
class IGAValues
{
 private:
	 /**
	 *pointer to NURBSMesh
	 */
	 NURBSMesh<dim>* mesh;
	 
	 /**
	 *quadrature point coordinates in real domain
	 */
   std::vector<std::vector<double> > QuadPoints_coords;
	 
	 /**
	 *quadrature point stencil (up to 9 points formula)
	 */
   std::vector<double> quadPtStencil;
	 
	 /**
	 *quadrature points rule used
	 */
   std::vector<std::vector<double> > defaultQuadPoints;
	 
	 /**
	 *basis function value at each quadrature point
	 */
   std::vector<std::vector<double> > basisValue;
	 
	 /**
	 *basis function first derivative at each quadrature point
	 */
   std::vector<std::vector<std::vector<double> > > basisDx;
	 
	 /**
	 *basis function second derivative at each quadrature point
	 */
   std::vector<std::vector<std::vector<std::vector<double> > > > basisDDx;
	 
	 /**
	 *determinant of jacobian times weight at each quadrature point (include quadpoints on plane and line) 
	 */
   std::vector<double> jacobian;
 public:
	 
	 /**
	 *IGAvalue constructor and destructor
	 */
   IGAValues(NURBSMesh<dim>* _mesh, unsigned int _dofPerControlPoint=0, unsigned int _numberOfDerivatives=0);
	 
	 /**
	 *update IGAValue for each cell
	 */
 	 void reinit(knotSpan<dim>& _cell, std::vector<std::vector<double> >* specialQuadPoints=0);
 
	 /**
	 *return coords of quadPoints in physical space
	 */
   std::vector<double>& quadraturePoints_coords( unsigned int q);
	
   /**
	 *return Jacobian times weight of volume
	 */
   double JxW(unsigned int q);
	 
   /**
	 *return Jacobian times weight of face
	 */
   double JxW_face(unsigned int q, unsigned int faceID);
	 
   /**
	 *return Jacobian times weight of line
	 */
   double JxW_line(unsigned int q, unsigned int lineID);
	 
   /**
	 *return Jacobian times weight of plane
	 */
   double JxW_plane(unsigned int q, unsigned int planeID);
	
	 /**
   *return basis function value at quadrature points of volume
	 */
   double shape_value(unsigned int i, unsigned int q);
	 
	 /**
   *return basis function value at quadrature points of face
	 */
   double shape_value_face(unsigned int i, unsigned int q, unsigned int faceID);
	 
	 /**
   *return basis function value at quadrature points of line
	 */
   double shape_value_line(unsigned int i, unsigned int q, unsigned int lineID);
	 
	 /**
   *return basis function value at quadrature points of plane
	 */
   double shape_value_plane(unsigned int i, unsigned int q, unsigned int planeID);
	
   /**
   *return basis function value first derivatives at quadrature points of volume
	 */
   std::vector<double>& shape_grad(unsigned int i, unsigned int q);
	 
   /**
   *return basis function value first derivatives at quadrature points of face
	 */
   std::vector<double>& shape_grad_face(unsigned int i, unsigned int q, unsigned int faceID);
	 
   /**
   *return basis function value first derivatives at quadrature points of line
	 */
   std::vector<double>& shape_grad_line(unsigned int i, unsigned int q, unsigned int lineID);
	 
   /**
   *return basis function value first derivatives at quadrature points of plane
	 */
   std::vector<double>& shape_grad_plane(unsigned int i, unsigned int q, unsigned int planeID);
  
	 /**
   *return basis function value second derivatives at quadrature points of volume
	 */
   std::vector<std::vector<double> >& shape_grad_grad(unsigned int i, unsigned int q);
	 
	 /**
   *return basis function value second derivatives at quadrature points of face
	 */
   std::vector<std::vector<double> >& shape_grad_grad_face(unsigned int i, unsigned int q, unsigned int faceID);
	 
	 /**
   *return basis function value second derivatives at quadrature points of line
	 */
   std::vector<std::vector<double> >& shape_grad_grad_line(unsigned int i, unsigned int q, unsigned int lineID);
	 
	 /**
   *return basis function value second derivatives at quadrature points of plane
	 */
   std::vector<std::vector<double> >& shape_grad_grad_plane(unsigned int i, unsigned int q, unsigned int planeID);
		
	 /**
	 *if cell contains boundary or marked line and plane, get # of quadpoint at boundary face
	 */
   std::vector<std::vector<unsigned int> > quadratureMap;
	 
	 /**
	 *if cell contains boundary or marked line and plane, get # of quadpoint along line
	 */
   std::vector<std::vector<unsigned int> > LinequadratureMap;	
	 
	 /**
	 *if cell contains boundary or marked line and plane, get # of quadpoint on plane
	 */
   std::vector<std::vector<unsigned int> > PlanequadratureMap;
  
	 /**
	 *return # of dof at control point at this # of dofs_per_cell
	 */
   unsigned int system_to_component_index(unsigned int dof);
	 
	 /**
	 *return # of control point at this # of dofs_per_cell
	 */
   unsigned int system_to_controlpoint_index(unsigned int dof);
	
	 /**
	 *pointer to knotSpan
	 *ofer another way to access knotSpan by IGAValues
	 */
   knotSpan<dim>* cell;
	 
	 /**
	 *number of quadrature_points of volume
	 */
   unsigned int n_quadrature_points;
	 
	 /**
	 *number of quadrature_points of boundary surface
	 */
   unsigned int n_face_quadrature_points;
	 
	 /**
	 *number of quadrature_points of marked line
	 */
   unsigned int n_line_quadrature_points;
	 
	 /**
	 *number of quadrature_points of marked plane
	 */
   unsigned int n_plane_quadrature_points;
	 
	 /**
	 *dofs_per_cell: total number of dofs of this cell
	 */
   unsigned int dofs_per_cell;
	 
	 /**
	 *controlpoints_per_cell: total controlpoints of this cell
	 */
   unsigned int controlpoints_per_cell;
	 
	 /**
	 *dofPerControlPoint: dof per control point
	 */
   unsigned int dofPerControlPoint;
	 
	 /**
	 *numberOfDerivatives: number Of Derivatives of basis function
	 */
   unsigned int numberOfDerivatives;
	 
	 /**
	 *quad point including additional quadpoints (face, line and plane) positions in physical space
	 */
   std::vector<std::vector<double> > quadPointLocations;
};

#endif
