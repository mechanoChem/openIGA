/*
*essential functions of IGAbase
*
*/
#include "../../../include/IGAbase/NURBSClass.h"

/*
*fill controlPointVector and knotSpanVector
*
*/
template <int dim>
NURBSMesh<dim>::NURBSMesh(nativeNURBSStructure<dim>& mesh, unsigned int _dofPerControlPoint, unsigned int _quadPtStencilSize): \
	dofPerControlPoint(_dofPerControlPoint), quadPtStencilSize(_quadPtStencilSize){	
	//fill controlPointVector
	unsigned int dofID=0;
	for (unsigned int controlPointID=0; controlPointID<mesh.ctrlPtCoefs[0].size(); controlPointID++){
		std::vector<unsigned int> dofs(dim);
		for (unsigned int a=0; a<dofPerControlPoint; a++) dofs[a]=dofID++;
		//create controlPoint and push back into controlPointVector
		controlPointVector.push_back(*(new controlPoint<dim>(controlPointID, mesh, dofs)));
	}
	//fill knotSpanVector
	IGAValues<dim> fe_values(this); //used to compute knotSpan edgeCoords and center 
	unsigned int numKnots=1, knotSpanID=0;
	for (unsigned int a=0; a<dim; a++) numKnots*=mesh.knots[a].size();
	for (unsigned int tempID=0; tempID<numKnots; tempID++){
		std::vector<unsigned int> tempIndices(dim);
		if (dim==3){
			std::div_t divresult= std::div ((int) tempID, (int) mesh.knots[0].size()*mesh.knots[1].size()); tempIndices[2]=divresult.quot;
			divresult= std::div ((int) divresult.rem, (int) mesh.knots[0].size()); tempIndices[1]=divresult.quot; tempIndices[0]=divresult.rem;
		}
		else if (dim==2){
			std::div_t divresult= std::div ((int) tempID, (int) mesh.knots[0].size()); 
			tempIndices[1]=divresult.quot; tempIndices[0]=divresult.rem;
		}
		else tempIndices[0]=tempID;
		
		//check if span width is greater than zero in each direction 
		bool skipLoopLast=false, skipLoopRep=false;
		for (unsigned int a=0; a<dim; a++) {
			if (tempIndices[a]==(mesh.knots[a].size()-1))  skipLoopLast=true;
			else if (mesh.knots[a].at(tempIndices[a])==mesh.knots[a].at(tempIndices[a]+1))  skipLoopRep=true;
		}
		if (skipLoopLast || skipLoopRep) continue;
		
		//create knotSpan and push back into knotSpanVector
		knotSpanVector.push_back(*(new knotSpan<dim>(knotSpanID++, mesh, *this, tempIndices, fe_values)));
	}
	printf("IGA data structures generated\n");
}

template <int dim>
NURBSMesh<dim>::~NURBSMesh (){}

template class NURBSMesh<1>;
template class NURBSMesh<2>;
template class NURBSMesh<3>;
