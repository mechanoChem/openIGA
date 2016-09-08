#include "../../include/supplementary/solutionClass.h"
/*
*class for solution storage
*data type include scalar, vector and tensor
*/
template <int dim>
solutionClass<dim>::solutionClass(NURBSMesh<dim>& _mesh, dataLocation _datalocation, dataType _datatype, std::string _variableName): mesh(&_mesh), datalocation(_datalocation), datatype(_datatype), variableName(_variableName){
	numQuadPoints=std::pow(mesh->quadPtStencilSize,dim);
	if (datalocation==NODAL){
		if (datatype==SCALAR) {values.resize(mesh->controlPointVector.size(), 0.0); numVariablesPerPoint=1;}
		else if (datatype==VECTOR) {values.resize(mesh->controlPointVector.size()*dim, 0.0); numVariablesPerPoint=dim;}
		else if (datatype==TENSOR) {values.resize(mesh->controlPointVector.size()*dim*dim, 0.0); numVariablesPerPoint=dim*dim;}
		else {printf("unknown dataType\n"); exit(-1);}
	}
	else if (datalocation==QUADRATURE){
		if (datatype==SCALAR) {values.resize(mesh->knotSpanVector.size()*numQuadPoints, 0.0); numVariablesPerPoint=1;}
		else if (datatype==VECTOR) {values.resize(mesh->knotSpanVector.size()*dim*numQuadPoints, 0.0); numVariablesPerPoint=dim;}
		else if (datatype==TENSOR) {values.resize(mesh->knotSpanVector.size()*dim*dim*numQuadPoints, 0.0); numVariablesPerPoint=dim*dim;}
		else {printf("unknown dataType\n"); exit(-1);}
		if (mass_sparsity_pattern==0) {mass_sparsity_pattern=new sparsityPattern(); mass_sparsity_pattern->init(mesh, true);}
		if (massMatrix==0) createMassMatrix();
	}
	else {printf("unknown dataLocation\n"); exit(-1);}
}

/*
*create mass jacobian
*/
template <int dim>
void solutionClass<dim>::createMassMatrix(){
	if (massMatrix) {printf("massMatrix already initialized\n"); exit(-1);}
	massMatrix=new sparseMatrix; 
	massMatrix->reinit(*mass_sparsity_pattern); *massMatrix=0;
	IGAValues<dim> fe_values(mesh, 1); 
		
	for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
		fe_values.reinit (*cell);
		unsigned int n_q_points= fe_values.n_quadrature_points;
		unsigned int dofs_per_cell=cell->local_mass_dof_indices.size();
		denseMatrix local_matrix(dofs_per_cell, dofs_per_cell);
			
		//Mass matrix
		for (unsigned int i=0; i<dofs_per_cell; ++i) {
			for (unsigned int j=0; j<dofs_per_cell; ++j){
				local_matrix(i,j)=0.0;
				for (unsigned int q=0; q<n_q_points; q++){
					local_matrix(i,j)+= fe_values.shape_value(i,q)*fe_values.shape_value(j,q)*fe_values.JxW(q); // Mij= Ni*Nj
				}
				(*massMatrix)(cell->local_mass_dof_indices[i], cell->local_mass_dof_indices[j])+=local_matrix(i,j); //Global Assembly
			}
		}
	}
	//printf("mass matrix generated\n");
}

/*
*project result at quadrature point to nodes using quarature points interpolation
*/
template <int dim>
void solutionClass<dim>::projectQuadratureValues(){

	if (datalocation!=QUADRATURE) {printf("datalocation!=QUADRATURE, so projection invalid\n"); exit(-1);}
	
	if (projectedValues.size()>0) projectedValues.clear();
	if (datatype==SCALAR) projectedValues.resize(mesh->controlPointVector.size(), 0.0);
	else if (datatype==VECTOR) projectedValues.resize(mesh->controlPointVector.size()*dim, 0.0);
	else if (datatype==TENSOR) projectedValues.resize(mesh->controlPointVector.size()*dim*dim, 0.0);
	
	//Create RHS
	denseVector system_rhs; system_rhs.reinit(*mass_sparsity_pattern);
	IGAValues<dim> fe_values(mesh, 1);
	//loop over all variables
	std::vector<double> tempValues(mesh->controlPointVector.size(), 0.0);
	for (unsigned int var=0; var<numVariablesPerPoint; var++){
		 system_rhs=0;
		 for (typename std::vector<knotSpan<dim> >::iterator cell=mesh->knotSpanVector.begin(); cell<mesh->knotSpanVector.end(); cell++){
			fe_values.reinit (*cell);
			unsigned int n_q_points= fe_values.n_quadrature_points;
			unsigned int dofs_per_cell=cell->local_mass_dof_indices.size();
			denseVector local_rhs(dofs_per_cell);
				
			//Mass matrix
			for (unsigned int i=0; i<dofs_per_cell; ++i) {
				local_rhs(i)=0.0;
				for (unsigned int q=0; q<n_q_points; q++){
					if (datatype==SCALAR) local_rhs(i)+= fe_values.shape_value(i,q)*((*this)(cell->id, q))*fe_values.JxW(q);
					else local_rhs(i)+= fe_values.shape_value(i,q)*((*this)(cell->id, q, var))*fe_values.JxW(q);
				}
				//if(var==0) std::cout<<"((*this)(cell->id, q, var)"<<(*this)(cell->id, 0, var)<<"fe_values.JxW(q)"<<fe_values.JxW(0)<<std::endl;
				//if(var==1) std::cout<<"eee="<<fe_values.shape_value(i,0)*((*this)(cell->id, 0, var))*fe_values.JxW(0)<<std::endl;
				system_rhs(cell->local_mass_dof_indices[i]) += local_rhs(i); //Global Assembly
			}
		} 
		//project
		for(std::vector<double>::iterator it=tempValues.begin(); it<tempValues.end(); it++) *it=0.0; //initializing solution vector
		int status=luSolver(&massMatrix->nonZeroValues.at(0), &mass_sparsity_pattern->columnIndices.at(0), &mass_sparsity_pattern->rowIndex.at(0), &system_rhs.values.at(0), (int) system_rhs.size(), (int) system_rhs.size(), (int) mass_sparsity_pattern->nnz, &tempValues.at(0), 0);
		if (status!=0) {printf("Solver exit status:%u\n", status); exit(-1);}
		for (unsigned int a=0; a<mesh->controlPointVector.size(); a++){
			projectedValues.at(a*numVariablesPerPoint+var)=tempValues.at(a);
		}
	}
}

template <int dim>
sparseMatrix* solutionClass<dim>::massMatrix=0;

template <int dim>
sparsityPattern* solutionClass<dim>::mass_sparsity_pattern=0;

template <int dim>
solutionClass<dim>::~solutionClass(){
	if (mass_sparsity_pattern) {delete mass_sparsity_pattern; mass_sparsity_pattern=0;}
	if (massMatrix) {delete massMatrix; massMatrix=0;}
}

template class solutionClass<1>;
template class solutionClass<2>;
template class solutionClass<3>;