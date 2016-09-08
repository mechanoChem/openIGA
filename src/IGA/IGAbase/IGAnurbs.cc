/*
*essential functions of IGAbase
*
*/

#include "../../../include/IGAbase/NURBSClass.h"

//basisFunction
template <int dim>
basisFunction<dim>::basisFunction(): weight(0.0) {p.resize(dim); localKnots.resize(dim);}


/*
*detect control point which belongs to the knotspan e.g. value is nonzero
*/
template <int dim>
bool basisFunction<dim>::pointInBasisSupport(std::vector<double>& coords){
	bool flag=true;
	for (unsigned int i=0; i<dim; i++){
		if  (!((coords[i]>=localKnots[i][0] && (coords[i]<=localKnots[i].back())))) flag=false;
	}
	return flag;
}

/*
*detect control point inside knot span
*/
template <int dim>
bool basisFunction<dim>::pointCrossSurfaceSupport(double x){
  bool flag=false;
  if  (x>localKnots[0][0] && x<localKnots[0].back()) flag=true;
  return flag;
}

/*
*recursive evaluation of value of NURBS using Cox-DeBoor recursion formula
*/
template <int dim>
double basisFunction<dim>::evaluateN(unsigned int i, unsigned int localKnotIndex, unsigned int _p, double xi, unsigned int k){
	double val=0.0;
	//Cox-DeBoor recursion formula
	if (_p==0){
		if ((xi>=localKnots[i][localKnotIndex]) && (xi<=localKnots[i][localKnotIndex+1]) && (k==0)) val=1;
		else val=0;
	}
	else{
		if (k>0){ 
			//recursion formula for derivative of basis function
			if (localKnots[i][localKnotIndex+_p]!=localKnots[i][localKnotIndex])     val+=   evaluateN(i, localKnotIndex, _p-1, xi, k-1)*(_p)/(localKnots[i][localKnotIndex+_p]-localKnots[i][localKnotIndex]);
			if (localKnots[i][localKnotIndex+_p+1]!=localKnots[i][localKnotIndex+1]) val+=-evaluateN(i, localKnotIndex+1, _p-1, xi, k-1)*(_p)/(localKnots[i][localKnotIndex+_p+1]-localKnots[i][localKnotIndex+1]);
		}
		else{
			//recursion formula for evaluation of basis function
			if (localKnots[i][localKnotIndex+_p]!=localKnots[i][localKnotIndex]) val+=evaluateN(i, localKnotIndex, _p-1, xi)*(xi-localKnots[i][localKnotIndex])/(localKnots[i][localKnotIndex+_p]-localKnots[i][localKnotIndex]);
			if (localKnots[i][localKnotIndex+_p+1]!=localKnots[i][localKnotIndex+1]) val+=evaluateN(i, localKnotIndex+1, _p-1, xi)*(localKnots[i][localKnotIndex+_p+1]-xi)/(localKnots[i][localKnotIndex+_p+1]-localKnots[i][localKnotIndex+1]);
		}
	}
	return val;
}

/*
*evaluate value of basisfunction:N*weight
*/
template <int dim>
double basisFunction<dim>::value(std::vector<double>& coords, std::vector<unsigned int>& dervs){ 
	double temp=1.0;
	for (unsigned int a=0; a<dim; a++){
		temp*=evaluateN(a, 0, p[a], coords[a], dervs[a]);
	}
	return temp*weight;
}

template class basisFunction<1>;
template class basisFunction<2>;
template class basisFunction<3>;

//controlPoint
/*
*process controal points from mesh
*/
template <int dim>
controlPoint<dim>::controlPoint(unsigned int _id, nativeNURBSStructure<dim>& mesh, std::vector<unsigned int>& _dofs): id(_id), dofs(_dofs), basis(){
	coords.resize(dim); dofs.resize(dim);
	//fill basis variables
	basis.weight=mesh.ctrlPtCoefs[dim].at(id);
	for (unsigned int a=0; a<dim; a++) {
		//fill control point variables
		coords[a]=mesh.ctrlPtCoefs[a].at(id)/basis.weight; //the division by basis.weight is required as octave nurbs outputs coords multiplied by the weights. Not required if nurbs coefficients obtained from any other source.
		//fill basis variables
		basis.p[a]=mesh.order[a];
		for (unsigned int b=0; b<mesh.order[a]+2; b++) basis.localKnots[a].push_back(mesh.knots[a].at(mesh.ctrlPtKnots[a][id]+b));
	}
}

template class controlPoint<1>;
template class controlPoint<2>;
template class controlPoint<3>;

/*
*process knot span from mesh
*storage contral points inside knot span
*generate local_dof_indices to global_dof_indices
*generate edge and center coords of knot span
*/
template <int dim>
knotSpan<dim>::knotSpan(unsigned int _id,  nativeNURBSStructure<dim>& mesh, NURBSMesh<dim>& nurbsmesh, std::vector<unsigned int>& tempIndices, IGAValues<dim>& fe_values): \
id(_id), boundaryFlags(2*dim, 0), lineFlags((dim*2*(dim-1))+1,0),planeFlags(2*dim+1, 0),defectFlags(6,0), Sp_lineQuad_Point(dim-1,0), endKnots(dim) {
	//fill endKnots
	for (unsigned int a=0; a<dim; a++) {endKnots[a].push_back(mesh.knots[a].at(tempIndices[a])); endKnots[a].push_back(mesh.knots[a].at(tempIndices[a]+1));}
	//fill controlPoints
	for (typename std::vector<controlPoint<dim> >::iterator a=nurbsmesh.controlPointVector.begin(); a<nurbsmesh.controlPointVector.end(); a++){ 
		std::vector<double> tempVec(dim);
		for (unsigned int b=0; b<dim; b++) tempVec[b]=0.5*(endKnots[b][0]+endKnots[b][1]);
		if (a->basis.pointInBasisSupport(tempVec)) {
			controlPoints.push_back(&(*a));
		      	if(a->basis.pointCrossSurfaceSupport(0.5)){
			  controlPoint_cross_surface_cof.push_back(0.5);
		      	}
		      	else controlPoint_cross_surface_cof.push_back(1.0);
		}
	}

	//fill local_dof_indices
	for (typename std::vector<controlPoint<dim>*>::iterator a=controlPoints.begin(); a<controlPoints.end(); a++){ 
		//fill local_dof_indices for this knotSpan
		for (unsigned int c=0; c<nurbsmesh.dofPerControlPoint; c++){
			unsigned int aIndex=((*a)->id)*nurbsmesh.dofPerControlPoint+c;
			local_dof_indices.push_back(aIndex);
		}
		local_mass_dof_indices.push_back((*a)->id);
	}
	//evaluate edgeCoords and center
	std::vector<double> quadPtStencil; 
	quadPtStencil.push_back(-1.0); quadPtStencil.push_back(1.0); quadPtStencil.push_back(1.0); quadPtStencil.push_back(1.0);
	//fill QuadPoints
	std::vector<std::vector<double> > quadPoints(std::pow(2, dim));
	std::vector<unsigned int> indices;
	for (unsigned int quadPointIndex=0; quadPointIndex<quadPoints.size(); quadPointIndex++){
		resolve<dim>(quadPointIndex, 2, indices); double temp=1.0;
		for (unsigned int j=0; j<dim; j++){
			quadPoints[quadPointIndex].push_back(quadPtStencil.at(2*indices[j]));
			temp*=quadPtStencil.at(2*indices[j]+1); 
		}
		quadPoints[quadPointIndex].push_back(temp);
	}
	fe_values.reinit(*this, &quadPoints);
	//fill edgeCoords and center
	edgeCoords.resize(quadPoints.size());
	center.resize(dim, 0.0);
	for (unsigned int quadPointIndex=0; quadPointIndex<quadPoints.size(); quadPointIndex++){
		for (unsigned int i=0; i<dim; i++) {
			edgeCoords[quadPointIndex].push_back(fe_values.quadPointLocations[quadPointIndex][i]);
			center[i]+=fe_values.quadPointLocations[quadPointIndex][i]/quadPoints.size();
		}
	}
}

/*
*determine wether points inside the knot span
*/
template <int dim>
bool knotSpan<dim>::pointInKnotSpan(std::vector<double>& coords){
	bool temp=true;
	for (unsigned int i=0; i<dim; i++){
		if  (!((coords[i]>=endKnots[i][0] && (coords[i]<=endKnots[i][1])))) temp=false;
	}
	return temp;
}

template class knotSpan<1>;
template class knotSpan<2>;
template class knotSpan<3>;

