#include "../../include/supplementary/sparsityDataClass.h"

sparsityPattern::sparsityPattern(): nnz(0), isMassMatrixPattern(false){}

void sparsityPattern::intialize(std::vector<knotSpan<1> >& knotSpanVector){
	//fill DOFConnections
	for (std::vector<knotSpan<1> >::iterator a=knotSpanVector.begin(); a<knotSpanVector.end(); a++){
		std::vector<unsigned int> local_dof_indices;
		if (isMassMatrixPattern) local_dof_indices=a->local_mass_dof_indices;
		else 					 local_dof_indices=a->local_dof_indices;
		for (std::vector<unsigned int>::iterator b=local_dof_indices.begin(); b<local_dof_indices.end(); b++){
			for (std::vector<unsigned int>::iterator c=local_dof_indices.begin(); c<local_dof_indices.end(); c++){
				DOFConnections[*b].insert(*c);
			}
		}
	}
	
	//fill nzMap, columnIndices, rowIndex
	rowIndex.resize(DOFConnections.size()+1, 0.0);
	nzMap.resize(DOFConnections.size());
	nnz=0; unsigned int DOFIndex=0;
	for (std::map<unsigned int, std::set<unsigned int> >::iterator a=DOFConnections.begin(); a!=DOFConnections.end(); a++) {
		rowIndex.at(DOFIndex)=nnz; //row index
		for (std::set<unsigned int>::iterator b=a->second.begin(); b!=a->second.end(); b++){
			columnIndices.push_back(*b); //column indices
			nzMap.at(DOFIndex)[*b]=nnz; //nzMap
			nnz++;
		}
		DOFIndex++;
	}
	rowIndex.back()=nnz; //this last entry is required for the SUPERLU NRFormat	
}

void sparsityPattern::intialize(std::vector<knotSpan<2> >& knotSpanVector){
	//fill DOFConnections
	for (std::vector<knotSpan<2> >::iterator a=knotSpanVector.begin(); a<knotSpanVector.end(); a++){
		std::vector<unsigned int> local_dof_indices;
		if (isMassMatrixPattern) local_dof_indices=a->local_mass_dof_indices;
		else local_dof_indices=a->local_dof_indices;
		for (std::vector<unsigned int>::iterator b=local_dof_indices.begin(); b<local_dof_indices.end(); b++){
			for (std::vector<unsigned int>::iterator c=local_dof_indices.begin(); c<local_dof_indices.end(); c++){
				DOFConnections[*b].insert(*c);
			}
		}
	}
	
	//fill nzMap, columnIndices, rowIndex
	rowIndex.resize(DOFConnections.size()+1, 0.0);
	nzMap.resize(DOFConnections.size());
	nnz=0; unsigned int DOFIndex=0;
	for (std::map<unsigned int, std::set<unsigned int> >::iterator a=DOFConnections.begin(); a!=DOFConnections.end(); a++) {
		rowIndex.at(DOFIndex)=nnz; //row index
		for (std::set<unsigned int>::iterator b=a->second.begin(); b!=a->second.end(); b++){
			columnIndices.push_back(*b); //column indices
			nzMap.at(DOFIndex)[*b]=nnz; //nzMap
			nnz++;
		}
		DOFIndex++;
	}
	rowIndex.back()=nnz; //this last entry is required for the SUPERLU NRFormat	
}

void sparsityPattern::intialize(std::vector<knotSpan<3> >& knotSpanVector){
	//fill DOFConnections
	for (std::vector<knotSpan<3> >::iterator a=knotSpanVector.begin(); a<knotSpanVector.end(); a++){
		std::vector<unsigned int> local_dof_indices;
		if (isMassMatrixPattern) local_dof_indices=a->local_mass_dof_indices;
		else local_dof_indices=a->local_dof_indices;
		for (std::vector<unsigned int>::iterator b=local_dof_indices.begin(); b<local_dof_indices.end(); b++){
			for (std::vector<unsigned int>::iterator c=local_dof_indices.begin(); c<local_dof_indices.end(); c++){
				DOFConnections[*b].insert(*c);
			}
		}
	}
	
	//fill nzMap, columnIndices, rowIndex
	rowIndex.resize(DOFConnections.size()+1, 0.0);
	nzMap.resize(DOFConnections.size());
	nnz=0; unsigned int DOFIndex=0;
	for (std::map<unsigned int, std::set<unsigned int> >::iterator a=DOFConnections.begin(); a!=DOFConnections.end(); a++) {
		rowIndex.at(DOFIndex)=nnz; //row index
		for (std::set<unsigned int>::iterator b=a->second.begin(); b!=a->second.end(); b++){
			columnIndices.push_back(*b); //column indices
			nzMap.at(DOFIndex)[*b]=nnz; //nzMap
			nnz++;
		}
		DOFIndex++;
	}
	rowIndex.back()=nnz; //this last entry is required for the SUPERLU NRFormat	
}

void sparsityPattern::init(NURBSMesh<1>* mesh, bool _isMassMatrixPattern){
	isMassMatrixPattern=_isMassMatrixPattern; intialize(mesh->knotSpanVector);
}
void sparsityPattern::init(NURBSMesh<2>* mesh, bool _isMassMatrixPattern){
	isMassMatrixPattern=_isMassMatrixPattern; intialize(mesh->knotSpanVector);
}
void sparsityPattern::init(NURBSMesh<3>* mesh, bool _isMassMatrixPattern){
	isMassMatrixPattern=_isMassMatrixPattern; intialize(mesh->knotSpanVector);
}

//*******
//*******

//denseMatrix
denseMatrix::denseMatrix(unsigned int _size1, unsigned int _size2){
	values.resize(_size1);
	for(std::vector<std::vector<double> >::iterator it=values.begin(); it<values.end(); it++) it->resize(_size2, 0.0);
}
unsigned int denseMatrix::size(){
	return values.size();
}
//operator  overloading
double& denseMatrix::operator() (unsigned int _a, unsigned int _b) {
	return values.at(_a).at(_b);
}
void denseMatrix::operator= (double _value) {
	if (_value==0.0){
		for(std::vector<std::vector<double> >::iterator it=values.begin(); it<values.end(); it++) {
			for (std::vector<double>::iterator b=it->begin(); b<it->end(); b++) *b=0.0;
		}
	}
	else{
		printf("denseMatrix operator '=' only legal when set to zero\n"); exit(1);
	}
}
void denseMatrix::print(){
	printf("\n");
	for (unsigned int i=0; i<values.size(); i++){
		for (unsigned int j=0; j<values.size(); j++){
			double value=values.at(i).at(j);
			printf("%8.1e ", value);
		}
		printf("\n");
	}
}

//*******
//*******

//denseVector
denseVector::denseVector(){}
denseVector::denseVector(unsigned int _size){
	reinit(_size);
}
unsigned int denseVector::size(){
	return values.size();
}
double denseVector::l2_norm(){
	double temp=0.0;
	for (std::vector<double>::iterator a=values.begin(); a<values.end(); a++) temp+=std::pow(*a, 2.0);
	return std::pow(temp, 0.5);
}
void denseVector::reinit(sparsityPattern& sparsitypattern){
	values.resize(sparsitypattern.nzMap.size(), 0.0);
}
void denseVector::reinit(unsigned int _size){
	values.resize(_size, 0.0);
}
//operator  overloading
double& denseVector::operator() (unsigned int _a) {
	return values.at(_a);}
void denseVector::operator= (double _value) {
	if (_value==0.0){
		for (std::vector<double>::iterator a=values.begin(); a<values.end(); a++) *a=0.0;
	}
	else{
		printf("denseVector operator '=' only legal when set to zero\n"); exit(1);
	}
}
void denseVector::operator= (denseVector& _value) {
	if (_value.size()==this->size()){
		unsigned int index=0;
		for (std::vector<double>::iterator a=values.begin(); a<values.end(); a++) *a=_value(index++);
	}
	else{
		printf("denseVector operator '=' only legal when both vectors of same size\n"); exit(1);
	}
}
void denseVector::operator+= (denseVector& _value) {
	if (_value.size()==this->size()){
		unsigned int index=0;
		for (std::vector<double>::iterator a=values.begin(); a<values.end(); a++) *a+=_value(index++);
	}
	else{
		printf("denseVector operator '+' only legal when both vectors of same size\n"); exit(1);
	}
}
void denseVector::print(){
	printf("\n");
	for (unsigned int i=0; i<this->size(); i++){
		printf("%8.1e ", values.at(i));
	}
	printf("\n");
}

//*******
//*******

//sparseMatrix
sparseMatrix::sparseMatrix(){sparsitypattern=0;}
void sparseMatrix::reinit(sparsityPattern& _sparsitypattern) {
	sparsitypattern=&_sparsitypattern; nonZeroValues.resize(sparsitypattern->nnz, 0.0);
}
//operator  overloading
double& sparseMatrix::operator() (unsigned int _a, unsigned int _b) {
	if (sparsitypattern->nzMap.at(_a).count(_b)>0) return nonZeroValues.at(sparsitypattern->nzMap.at(_a)[_b]);
	else {printf("nzMap error at a=%u, b=%u\n",_a,_b); exit(1);}
}
double sparseMatrix::val (unsigned int _a, unsigned int _b) {
	if (sparsitypattern->nzMap.at(_a).count(_b)>0) return nonZeroValues.at(sparsitypattern->nzMap.at(_a)[_b]);
	else { return 0.0;}
}
void sparseMatrix::operator= (double _value) {
	if (_value==0.0){
		for (std::vector<double>::iterator a=nonZeroValues.begin(); a<nonZeroValues.end(); a++) *a=0.0;
	}
	else{
		printf("sparseMatrix operator '=' only legal when set to zero\n"); exit(1);
	}
}
void sparseMatrix::print(){
	printf("\n");
	for (unsigned int i=0; i<sparsitypattern->nzMap.size(); i++){
		for (unsigned int j=0; j<sparsitypattern->nzMap.size(); j++){
			double value=0.0;
			if (sparsitypattern->nzMap.at(i).count(j)>0) value=(*this)(i,j);
			printf("%8.1e ", value);
		}
		printf("\n");
	}
}
