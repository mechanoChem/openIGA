/*
*essential functions of IGAbase
*/
#include "../../../include/IGAbase/NURBSClass.h"
#include <H5Cpp.h>

using namespace std;
/*
*Initialize class nativeNURBSStructure
*
*/
template <int dim>
nativeNURBSStructure<dim>::nativeNURBSStructure(const char* _hdf5FileName){
	order.resize(dim);
	numCtrlPts.resize(dim);
	knots.resize(dim);
	ctrlPtCoefs.resize(dim+1);//last one store the weight of control point
	ctrlPtKnots.resize(dim);
	hdf5FileName=_hdf5FileName;
	readHDF5();
}

template <int dim>
nativeNURBSStructure<dim>::~nativeNURBSStructure (){}

/*
*read .h5 mesh and fill mesh information
*
*/
template <int dim>
int nativeNURBSStructure<dim>::readHDF5(){
	if (dim<0 || dim>3){printf("Error: (dim<0 || dim>3). dim=%u\n", dim); exit(-1);}
	
	//Opening NURBS mesh file
	H5::H5File file (hdf5FileName, H5F_ACC_RDONLY);
	
	//Opening intial groups to access the NURBS datasets
	H5::Group grp_mesh  = file.openGroup("mesh"); 
	H5::Group grp_value  = grp_mesh.openGroup("value");
	
	//polynomial order along each direction
	H5::Group grp_order  = grp_value.openGroup("order");
	H5::DataSet dataset_order = grp_order.openDataSet("value");
	H5::DataSpace dataspace_order = dataset_order.getSpace();
	int rank = dataspace_order.getSimpleExtentNdims();
	hsize_t dims_order[2]; dataspace_order.getSimpleExtentDims(dims_order, NULL);
	if ( (rank>0 && dims_order[0]!=dim) || (rank==0 && dim!=1) ){printf("Error: Mesh object and mesh file have different dimensions. HDF5 file dim:%u, but specified dim:%u\n", (int) dims_order[0], dim); exit(-1);}
	unsigned int data_order[dim];
	dataset_order.read(data_order, H5::PredType::NATIVE_INT);
	for(unsigned int i=0; i<dim; i++){order[i]=data_order[i]-1;}
	//printf("Size is:%u\n%u %u %u\n", mesh.order.size(),  mesh.order.at(0),  mesh.order.at(1),  mesh.order.at(2));
	grp_order.close();
	
	//number of control points along each direction
	H5::Group grp_number  = grp_value.openGroup("number");
	H5::DataSet dataset_number = grp_number.openDataSet("value");
	unsigned int data_number[dim];
	dataset_number.read(data_number, H5::PredType::NATIVE_INT);
	for(unsigned int i=0; i<dim; i++){numCtrlPts[i]=data_number[i];}
	//printf("Size is:%u\n%u %u %u\n", mesh.numCtrlPts.size(),  mesh.numCtrlPts.at(0),  mesh.numCtrlPts.at(1),  mesh.numCtrlPts.at(2));
	grp_number.close();
			
	//knot values along each direction
	H5::Group grp_knots  = grp_value.openGroup("knots");
	if (dim==1){
		H5::DataSet dataset_knots_value = grp_knots.openDataSet("value");
		unsigned int num_knots=order.at(0)+numCtrlPts.at(0)+1; //number of knots = n+p+1
		double* data_knots= new double[num_knots];
		dataset_knots_value.read(data_knots, H5::PredType::NATIVE_DOUBLE);
		knots.at(0).resize(num_knots);
		for(unsigned int j=0; j<num_knots; j++){knots.at(0).at(j)=data_knots[j];}
		delete [] data_knots;
		dataset_knots_value.close();
	}
	else{
		H5::Group grp_knots_value = grp_knots.openGroup("value");
		for(unsigned int i=0; i<dim; i++){
			char buffer[10]; sprintf(buffer, "_%u", i);
			H5::Group grp_knots_value_i = grp_knots_value.openGroup(buffer);
			H5::DataSet dataset_knots_i = grp_knots_value_i.openDataSet("value");
			unsigned int num_knots=order.at(i)+numCtrlPts.at(i)+1; //number of knots = n+p+1
			double* data_knots= new double[num_knots]; 
			dataset_knots_i.read(data_knots, H5::PredType::NATIVE_DOUBLE);
			knots.at(i).resize(num_knots);
			for(unsigned int j=0; j<num_knots; j++){knots.at(i).at(j)=data_knots[j];}
			delete [] data_knots;
			grp_knots_value_i.close();
		}
	}
	grp_knots.close();
	
	//control point coordinates along each direction and their weights
	H5::Group grp_ctrlpts  = grp_value.openGroup("coefs");
	H5::DataSet dataset_ctrlpts = grp_ctrlpts.openDataSet("value");
	H5::DataSpace dataspace_ctrlpts = dataset_ctrlpts.getSpace();
	hsize_t dims_ctrlpts[dim+1]; dataspace_ctrlpts.getSimpleExtentDims(dims_ctrlpts, NULL);
	/*
	for(unsigned int i=0; i<dim+1; i++){
		printf("%u  ", dims_ctrlpts[i]);
	}
	printf("\n");
	*/
	//read control point coefficients
	unsigned int numPoints=1, numValues=1;
	for(unsigned int i=0; i<dim; i++){numPoints*=dims_ctrlpts[i];}
	for(unsigned int i=0; i<dim+1; i++){numValues*=dims_ctrlpts[i];}
	double* data_ctrlpts= new double [numValues];
	dataset_ctrlpts.read(data_ctrlpts, H5::PredType::NATIVE_DOUBLE);
	if (dim==3){
		//populate the ctrlPtCoefs vector
		for(unsigned int i=0; i<dims_ctrlpts[3]; i++){
			for(unsigned int j=0; j<dims_ctrlpts[0]; j++){
				for(unsigned int k=0; k<dims_ctrlpts[1]; k++){
					for(unsigned int l=0; l<dims_ctrlpts[2]; l++){
						ctrlPtCoefs.at(i).push_back(data_ctrlpts[j*dims_ctrlpts[1]*dims_ctrlpts[2]*dims_ctrlpts[3] + k*dims_ctrlpts[2]*dims_ctrlpts[3] + l*dims_ctrlpts[3] + i]);
					}
				}
			}
		}
		//populate the ctrlPtKnots vector
		for(unsigned int j=0; j<dims_ctrlpts[0]; j++){
			for(unsigned int k=0; k<dims_ctrlpts[1]; k++){
				for(unsigned int l=0; l<dims_ctrlpts[2]; l++){
					ctrlPtKnots.at(0).push_back(l); //x
					ctrlPtKnots.at(1).push_back(k); //y
					ctrlPtKnots.at(2).push_back(j); //z
				}
			}
		}
	}
	else if (dim==2){
		//populate the ctrlPtCoefs vector
		for(unsigned int i=0; i<dims_ctrlpts[2]; i++){
			for(unsigned int j=0; j<dims_ctrlpts[0]; j++){
				for(unsigned int k=0; k<dims_ctrlpts[1]; k++){
					if (i==2) continue; //skip z coodinates
					if (i==3) ctrlPtCoefs.at(i-1).push_back(data_ctrlpts[j*dims_ctrlpts[1]*dims_ctrlpts[2] + k*dims_ctrlpts[2] + i]);
					else ctrlPtCoefs.at(i).push_back(data_ctrlpts[j*dims_ctrlpts[1]*dims_ctrlpts[2] + k*dims_ctrlpts[2] + i]);
				}
			}
		}
		//populate the ctrlPtKnots vector
		for(unsigned int j=0; j<dims_ctrlpts[0]; j++){
			for(unsigned int k=0; k<dims_ctrlpts[1]; k++){
				ctrlPtKnots.at(0).push_back(k); //x
				ctrlPtKnots.at(1).push_back(j); //y
			}
		}
	}
	else{
		//populate the ctrlPtCoefs vector
		for(unsigned int i=0; i<dims_ctrlpts[1]; i++){
			for(unsigned int j=0; j<dims_ctrlpts[0]; j++){
				if (i==1 || i==2) continue; //skip y,z coodinates
				if (i==3) ctrlPtCoefs.at(i-2).push_back(data_ctrlpts[j*dims_ctrlpts[1] + i]);
				else ctrlPtCoefs.at(i).push_back(data_ctrlpts[j*dims_ctrlpts[1] + i]);
			}
		}
		//populate the ctrlPtKnots vector
		for(unsigned int j=0; j<dims_ctrlpts[0]; j++){
			ctrlPtKnots.at(0).push_back(j); //x
		}
	}
	delete [] data_ctrlpts;
	grp_ctrlpts.close();
	/*
	for(unsigned int i=0; i<dims_ctrlpts[3]; i++){
		for(unsigned int j=0; j<dims_ctrlpts[2]; j++){
			printf("%12.7f  ", data_ctrlpts[0][0][j][i]);
		}
		printf("\n");
	}
	printf("\n");
	*/
	
	/*for(unsigned int i=0; i<numPoints; i++){
		printf("%12.7f  ", mesh.ctrlPtCoefs.at(3).at(i));
	}
	printf("\n");
	*/
	file.close();
	printf("geometry read from:%s\n", hdf5FileName);
	return 1;
}

template class nativeNURBSStructure<1>;
template class nativeNURBSStructure<2>;
template class nativeNURBSStructure<3>;