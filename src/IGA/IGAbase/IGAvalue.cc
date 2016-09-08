/*
*essential functions of IGAbase
*
*/
#include "../../../include/IGAbase/NURBSClass.h"

using namespace std;

/*
*define quadrature rules
*
*/
template <int dim>
IGAValues<dim>::IGAValues(NURBSMesh<dim>* _mesh, unsigned int _dofPerControlPoint, unsigned int _numberOfDerivatives): \
mesh(_mesh), cell(0), n_quadrature_points(0), n_face_quadrature_points(0), dofs_per_cell(0), numberOfDerivatives(_numberOfDerivatives){
	if (_dofPerControlPoint==0) {dofPerControlPoint=mesh->dofPerControlPoint;} //default number of dofPerControlPoint
  else {dofPerControlPoint=_dofPerControlPoint;}
			
  //First derivative computation is always enabled, as it is used in the jacobian computation
  if (numberOfDerivatives==0) numberOfDerivatives=1;
  if (numberOfDerivatives>2) {printf("numberOfDerivatives>2 not yet implemented\n"); exit(-1);}
			
  //initialize quad points
  unsigned int quadPtStencilSize=mesh->quadPtStencilSize;
  if (quadPtStencilSize==1) {quadPtStencil.push_back(0.0); quadPtStencil.push_back(2.0);}
  else if (quadPtStencilSize==2) {quadPtStencil.push_back(-1.0/std::pow(3, 0.5)); quadPtStencil.push_back(1.0); \
    quadPtStencil.push_back( 1.0/std::pow(3, 0.5)); quadPtStencil.push_back(1.0);}
  else if (quadPtStencilSize==3) {quadPtStencil.push_back(0.0); quadPtStencil.push_back(8.0/9.0); \
    quadPtStencil.push_back(-std::pow(3.0/5.0, 0.5)); quadPtStencil.push_back(5.0/9.0); \
    quadPtStencil.push_back( std::pow(3.0/5.0, 0.5)); quadPtStencil.push_back(5.0/9.0);}
  else if (quadPtStencilSize==4) {quadPtStencil.push_back( std::pow((3.0-2.0*std::pow(6.0/5.0, 0.5))/7.0, 0.5)); quadPtStencil.push_back((18.0+std::pow(30.0,0.5))/36.0); \
    quadPtStencil.push_back(-std::pow((3.0-2.0*std::pow(6.0/5.0, 0.5))/7.0, 0.5)); quadPtStencil.push_back((18.0+std::pow(30.0,0.5))/36.0); \
    quadPtStencil.push_back( std::pow((3.0+2.0*std::pow(6.0/5.0, 0.5))/7.0, 0.5)); quadPtStencil.push_back((18.0-std::pow(30.0,0.5))/36.0); \
    quadPtStencil.push_back(-std::pow((3.0+2.0*std::pow(6.0/5.0, 0.5))/7.0, 0.5)); quadPtStencil.push_back((18.0-std::pow(30.0,0.5))/36.0);}
  else if (quadPtStencilSize==5) {quadPtStencil.push_back(0.0); quadPtStencil.push_back(128.0/225.0); \
    quadPtStencil.push_back( std::pow(5.0-2.0*std::pow(10.0/7.0, 0.5), 0.5)/3.0); quadPtStencil.push_back((322.0+13.0*std::pow(70.0,0.5))/900.0); \
    quadPtStencil.push_back(-std::pow(5.0-2.0*std::pow(10.0/7.0, 0.5), 0.5)/3.0); quadPtStencil.push_back((322.0+13.0*std::pow(70.0,0.5))/900.0); \
    quadPtStencil.push_back( std::pow(5.0+2.0*std::pow(10.0/7.0, 0.5), 0.5)/3.0); quadPtStencil.push_back((322.0-13.0*std::pow(70.0,0.5))/900.0); \
    quadPtStencil.push_back(-std::pow(5.0+2.0*std::pow(10.0/7.0, 0.5), 0.5)/3.0); quadPtStencil.push_back((322.0-13.0*std::pow(70.0,0.5))/900.0);}
  else if (quadPtStencilSize==6) {
    quadPtStencil.push_back(0.6612093864662645); quadPtStencil.push_back(0.3607615730481386);
    quadPtStencil.push_back(-0.6612093864662645); quadPtStencil.push_back(0.3607615730481386);
    quadPtStencil.push_back(0.2386191860831969); quadPtStencil.push_back(0.4679139345726910);
    quadPtStencil.push_back(-0.2386191860831969); quadPtStencil.push_back(0.4679139345726910);
    quadPtStencil.push_back(0.9324695142031521); quadPtStencil.push_back(0.1713244923791704);
    quadPtStencil.push_back(-0.9324695142031521); quadPtStencil.push_back(0.1713244923791704);
  }
  else if (quadPtStencilSize==7) {quadPtStencil.push_back(0.0); quadPtStencil.push_back(0.4179591836734694);
    quadPtStencil.push_back(0.4058451513773972); quadPtStencil.push_back(0.3818300505051189);
    quadPtStencil.push_back(-0.4058451513773972); quadPtStencil.push_back(0.3818300505051189);
    quadPtStencil.push_back(0.7415311855993945); quadPtStencil.push_back(0.2797053914892766);
    quadPtStencil.push_back(-0.7415311855993945); quadPtStencil.push_back(0.2797053914892766);
    quadPtStencil.push_back(0.9491079123427585); quadPtStencil.push_back(0.1294849661688697);
    quadPtStencil.push_back(-0.9491079123427585); quadPtStencil.push_back(0.1294849661688697);
  }
  else if (quadPtStencilSize==8) {
    quadPtStencil.push_back(0.1834346424956498); quadPtStencil.push_back(0.3626837833783620); 
    quadPtStencil.push_back(-0.1834346424956498); quadPtStencil.push_back(0.3626837833783620);
    quadPtStencil.push_back(0.5255324099163290); quadPtStencil.push_back(0.3137066458778873);
    quadPtStencil.push_back(-0.5255324099163290); quadPtStencil.push_back(0.3137066458778873);
    quadPtStencil.push_back(0.7966664774136267); quadPtStencil.push_back(0.2223810344533745);
    quadPtStencil.push_back(-0.7966664774136267); quadPtStencil.push_back(0.2223810344533745);
    quadPtStencil.push_back(0.9602898564975363); quadPtStencil.push_back(0.1012285362903763);
    quadPtStencil.push_back(-0.9602898564975363); quadPtStencil.push_back(0.1012285362903763);
	}
  else if (quadPtStencilSize==9) {quadPtStencil.push_back(0.0); quadPtStencil.push_back(0.3302393550012598);
    quadPtStencil.push_back(-0.8360311073266358); quadPtStencil.push_back(0.1806481606948574);
    quadPtStencil.push_back(0.8360311073266358); quadPtStencil.push_back(0.1806481606948574);
    quadPtStencil.push_back(-0.9681602395076261); quadPtStencil.push_back(0.0812743883615744);
    quadPtStencil.push_back(0.9681602395076261); quadPtStencil.push_back(0.0812743883615744);
    quadPtStencil.push_back(-0.3242534234038089); quadPtStencil.push_back(0.3123470770400029);
    quadPtStencil.push_back(0.3242534234038089); quadPtStencil.push_back(0.3123470770400029);
    quadPtStencil.push_back(-0.6133714327005904); quadPtStencil.push_back(0.2606106964029354);
    quadPtStencil.push_back(0.6133714327005904); quadPtStencil.push_back(0.2606106964029354);
  }
  else {printf("%u point gauss quadrature not yet implemented\n", quadPtStencilSize); exit(-1);}
			
  //fill defaultQuadPoints
  defaultQuadPoints.resize(std::pow(quadPtStencilSize, dim));
  std::vector<unsigned int> indices;
  for (unsigned int quadPointIndex=0; quadPointIndex<defaultQuadPoints.size(); quadPointIndex++){
    resolve<dim>(quadPointIndex, quadPtStencilSize, indices); double temp=1.0;
    for (unsigned int j=0; j<dim; j++){
		  defaultQuadPoints[quadPointIndex].push_back(quadPtStencil.at(2*indices[j]));
		  temp*=quadPtStencil.at(2*indices[j]+1); 
    }
    defaultQuadPoints[quadPointIndex].push_back(temp);
  }
}
/*
*update IGAValues for given knotspan
*fill quadrature points inside the knotspan
* quadrature points include regular quadpoints, quad points on boundary face, on plane inside cell and along line inside cell
*calculate values of basis functions and its first and second derivatives
*calculate Joccobian for volume, face area and line length
*/
template <int dim>
void IGAValues<dim>::reinit(knotSpan<dim>& _cell, std::vector<std::vector<double> >* specialQuadPoints){		
  //Set the quadrature information (either default or special)
  std::vector<std::vector<double> > quadPoints;
  quadratureMap.clear(); quadratureMap.resize(2*dim); //2*dim for each face quad points
  LinequadratureMap.clear(); LinequadratureMap.resize(dim*2*(dim-1)); //
  PlanequadratureMap.clear();PlanequadratureMap.resize(2*dim);
 
  n_face_quadrature_points=0;//for boundary surface
  n_line_quadrature_points=0;
  n_plane_quadrature_points=0;
		
  if (specialQuadPoints) quadPoints=*specialQuadPoints;
  else quadPoints=defaultQuadPoints;//=pow(quadPtStencilSize, dim)
  n_quadrature_points=quadPoints.size(); //includes only cell quad points.

  //fill face quad points, if any face boundaryFlags marked
  std::map<unsigned int, unsigned int> tempFaceQuadMap;//for boundaryface
  std::map<unsigned int, unsigned int> tempLineQuadMap;//for boundary line
  std::map<unsigned int, unsigned int> tempPlaneQuadMap;//for inside plane
  if (!specialQuadPoints){
    //fill face quad points, if any face boundaryFlags marked
    for (unsigned int faceID=0; faceID<2*dim; faceID++) {
			if (_cell.boundaryFlags[faceID]>0){
  			if (!n_face_quadrature_points) n_face_quadrature_points=std::pow(mesh->quadPtStencilSize, dim-1);
  			std::vector<unsigned int> indices;
  			for (unsigned int quadPointIndex=0; quadPointIndex<std::pow(mesh->quadPtStencilSize, dim-1); quadPointIndex++){
    			quadratureMap[faceID].push_back(quadPoints.size());//# of quadPoint will set
    			tempFaceQuadMap[quadPoints.size()]=faceID/2;//quadPoints.size() is this # of quadPoint will set. faceID/2 to get dim of this face.
    			std::vector<double> tempQuadpts;
    			resolve<dim-1>(quadPointIndex, mesh->quadPtStencilSize, indices); 
    			double temp=1.0; unsigned int tempCount=0;
    			for (unsigned int j=0; j<dim; j++){
      			if (j==(faceID/2)) {tempQuadpts.push_back(std::pow(-1,(faceID+1)%2)); }//quad_point for normal direction -1 or 1.  
      			else {tempQuadpts.push_back(quadPtStencil.at(2*indices[tempCount])); temp*=quadPtStencil.at(2*indices[tempCount]+1); tempCount++;} 
    			}
    			tempQuadpts.push_back(temp);//temp is weights of quad_point
    			quadPoints.push_back(tempQuadpts);
  			}
			}
    }
    //fill line quad points and line quad points, if any lineFlags marked
    for (unsigned int lineID=0; lineID<(dim*2*(dim-1)); lineID++) {
			if (_cell.lineFlags[lineID]>0){
  			if (!n_line_quadrature_points) n_line_quadrature_points=mesh->quadPtStencilSize;
	  		for (unsigned int quadPointIndex=0; quadPointIndex<n_line_quadrature_points; quadPointIndex++){
    			LinequadratureMap[lineID].push_back(quadPoints.size());//index of this line quadPoint
    			tempLineQuadMap[quadPoints.size()]=lineID/(2*(dim-1));//determine direction(x,y,z) of this line quadPoint
    			std::vector<double> tempLineQuadpts;
    			int TemK=0;
    			for (unsigned int j=0; j<dim; j++){
     			 if (j==(lineID/(2*(dim-1)))) {tempLineQuadpts.push_back(quadPtStencil.at(2*quadPointIndex));}
     				else if(_cell.lineFlags[dim*2*(dim-1)]==0){
							int temDimID=(lineID+1)%(2*(dim-1));
							if(j<1 or dim==2) tempLineQuadpts.push_back(std::pow(-1,temDimID));
							if(dim==3 and j>=1){//push the thirds coords if necessary
	  					if(temDimID==1 or temDimID==2) {tempLineQuadpts.push_back(-1);}
	  					if(temDimID==3 or temDimID==0){tempLineQuadpts.push_back(1); }
        			}
      			}
      			else if(_cell.lineFlags[dim*2*(dim-1)]==1){
							tempLineQuadpts.push_back(_cell.Sp_lineQuad_Point[TemK]);
							TemK++;
      			}
   		 		}
         	tempLineQuadpts.push_back(quadPtStencil.at(2*quadPointIndex+1));//push weights
         	quadPoints.push_back(tempLineQuadpts);//push this line quadpoint
  		 	}
      }
    }
    
    //fill plane quad points, if any plane planeFlags marked
    for (unsigned int planeID=0; planeID<2*dim; planeID++) {
			if (_cell.planeFlags[planeID]>0){
  			if (!n_plane_quadrature_points) n_plane_quadrature_points=std::pow(mesh->quadPtStencilSize, dim-1);
  			std::vector<unsigned int> indices;
  			for (unsigned int quadPointIndex=0; quadPointIndex<n_plane_quadrature_points; quadPointIndex++){
    			PlanequadratureMap[planeID].push_back(quadPoints.size());//# of quadPoint will set
    			tempPlaneQuadMap[quadPoints.size()]=planeID/2;//quadPoints.size() is this # of quadPoint will set. faceID/2 to get dim of this face.
    			std::vector<double> tempPlaneQuadpts;
    			resolve<dim-1>(quadPointIndex, mesh->quadPtStencilSize, indices); 
    			double temp=1.0; unsigned int tempCount=0;
    			for (unsigned int j=0; j<dim; j++){
      			if (j==(planeID/2)) {
        			if(_cell.planeFlags[2*dim]==0) {tempPlaneQuadpts.push_back(std::pow(-1,(planeID+1)%2));}
							else if(_cell.planeFlags[2*dim]==1) {tempPlaneQuadpts.push_back(_cell.Sp_planeQuad_Point);}
    				}
      			else {tempPlaneQuadpts.push_back (quadPtStencil.at(2*indices[tempCount])); temp*=quadPtStencil.at(2*indices[tempCount]+1); tempCount++;} 
    			}
    			tempPlaneQuadpts.push_back(temp);//temp is weights of quad_point
    			quadPoints.push_back(tempPlaneQuadpts);
  			}
			}
    }    
  }

  //cell level calculations
  unsigned int temp_n_quadrature_points=quadPoints.size(); //may include face quad points in addition to cell quad points
  cell=&_cell;
  controlpoints_per_cell=cell->controlPoints.size();
  dofs_per_cell=cell->controlPoints.size()*dofPerControlPoint; //not including special DOF like lagrange multiplier
  jacobian.resize(temp_n_quadrature_points, 0.0);
  quadPointLocations.clear(); quadPointLocations.resize(temp_n_quadrature_points);
  basisValue.clear(); basisValue.resize(temp_n_quadrature_points); 
  basisDx.clear(); basisDx.resize(temp_n_quadrature_points); 
  if (numberOfDerivatives==2) {basisDDx.clear(); basisDDx.resize(temp_n_quadrature_points);}
		
  //Compute basis gradient values
  //Set size of shapeValue and shapeGrad
  std::vector<std::vector<std::vector<double> > > R1D;
  std::vector<std::vector<std::vector<std::vector<double> > > > R2D;
  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > R3D;
  if (dim==1) R1D.resize(temp_n_quadrature_points); 
  if (dim==2) R2D.resize(temp_n_quadrature_points); 
  if (dim==3) R3D.resize(temp_n_quadrature_points);
  unsigned int quadIndex, basisIndex; double tempVal;
  QuadPoints_coords.clear(); QuadPoints_coords.resize(temp_n_quadrature_points);
  for (quadIndex=0; quadIndex<temp_n_quadrature_points; quadIndex++) {
    basisValue[quadIndex].resize(controlpoints_per_cell, 0.0);
    std::vector<double> coords(dim, 0.0);
    for (unsigned int i=0; i<dim; i++) {
			coords[i]=cell->endKnots[i][0] +  0.5*(1+quadPoints.at(quadIndex)[i])*(cell->endKnots[i][1]-cell->endKnots[i][0]);
			QuadPoints_coords[quadIndex].push_back(coords[i]);
    }				
    if (dim==1){
			//Compute A's and W's
			std::vector<std::vector<double> > A(controlpoints_per_cell);
			std::vector<double> W(numberOfDerivatives+1, 0.0); 
			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
				A[basisIndex].resize(numberOfDerivatives+1, 0.0); 
				for (unsigned int l=0; l<=numberOfDerivatives; l++) {
					std::vector<unsigned int> dervs(dim, 0); dervs[0]=l;
					tempVal=cell->controlPoints.at(basisIndex)->basis.value(coords, dervs);
					A[basisIndex][l]=tempVal; W[l]+=tempVal;
				}
			}
			//Compute R's
			R1D[quadIndex].resize(controlpoints_per_cell);
			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
				R1D[quadIndex][basisIndex].resize(numberOfDerivatives+1);
				for (unsigned int l=0; l<=numberOfDerivatives; l++) {
					double v = 0.0; 
					for (unsigned int i=0; i<=l; i++) {
						if (i==0) v=A[basisIndex][l];
						else v -= Bin(l,i)*W[i]*R1D[quadIndex][basisIndex][l-i];
					}
					R1D[quadIndex][basisIndex][l]=v/W[0];
				}
			}
    }
    else if (dim==2){
			//Compute A's and W's
			std::vector<std::vector<std::vector<double> > > A(controlpoints_per_cell);
			std::vector<std::vector<double> > W(numberOfDerivatives+1); 
			for (unsigned int l=0; l<=numberOfDerivatives; l++) {W[l].resize(numberOfDerivatives-l+1, 0.0);}
			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
				A[basisIndex].resize(numberOfDerivatives+1); 
				for (unsigned int l=0; l<=numberOfDerivatives; l++) {
					A[basisIndex][l].resize(numberOfDerivatives-l+1, 0.0);
					for (unsigned int m=0; m<=numberOfDerivatives-l; m++) {
 					 	std::vector<unsigned int> dervs(dim, 0); dervs[0]=l; dervs[1]=m;
 						tempVal=cell->controlPoints.at(basisIndex)->basis.value(coords, dervs);
 				 	 A[basisIndex][l][m]=tempVal; W[l][m]+=tempVal;
 					}
				}
			}
			//Compute R's
			R2D[quadIndex].resize(controlpoints_per_cell);
			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
				R2D[quadIndex][basisIndex].resize(numberOfDerivatives+1);
  			for (unsigned int l=0; l<=numberOfDerivatives; l++) {
    			R2D[quadIndex][basisIndex][l].resize(numberOfDerivatives-l+1, 0.0);
    			for (unsigned int m=0; m<=numberOfDerivatives-l; m++) {
      			double v = 0.0;
      			for (unsigned int i=0; i<=l; i++) {
							double v1 = 0.0;
							for (unsigned int j=0; j<=m; j++) {
	 						 if ((i+j)==0) v  = A[basisIndex][l][m];  //i=j=0
	  					 else		v1 += Bin(m,j)*W[i][j]*R2D[quadIndex][basisIndex][l-i][m-j];
						 	}
							v -= Bin(l,i)*v1;
      			}
      			R2D[quadIndex][basisIndex][l][m]=v/W[0][0];
    			}
  			}
			}
    }
    else if (dim==3){
			//Compute A's and W's
			std::vector<std::vector<std::vector<std::vector<double> > > >A(controlpoints_per_cell);
			std::vector<std::vector<std::vector<double> > > W(numberOfDerivatives+1); 
			for (unsigned int l=0; l<=numberOfDerivatives; l++){
				W[l].resize(numberOfDerivatives-l+1);
				for (unsigned int m=0; m<=numberOfDerivatives-l; m++) W[l][m].resize(numberOfDerivatives-l-m+1, 0.0);
			}	
			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
				A[basisIndex].resize(numberOfDerivatives+1); 
				for (unsigned int l=0; l<=numberOfDerivatives; l++) {
					A[basisIndex][l].resize(numberOfDerivatives-l+1);
    			for (unsigned int m=0; m<=numberOfDerivatives-l; m++) {
      			A[basisIndex][l][m].resize(numberOfDerivatives-l-m+1, 0.0);
      			for (unsigned int n=0; n<=numberOfDerivatives-l-m; n++) {
							std::vector<unsigned int> dervs(dim, 0); 
							dervs[0]=l; dervs[1]=m; dervs[2]=n;
							tempVal=cell->controlPoints.at(basisIndex)->basis.value(coords, dervs);
							A[basisIndex][l][m][n]=tempVal; W[l][m][n]+=tempVal;
      			}
    			}
  			}
			}
			//Compute R's
			R3D[quadIndex].resize(controlpoints_per_cell);
			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
  			R3D[quadIndex][basisIndex].resize(numberOfDerivatives+1);
  			for (unsigned int l=0; l<=numberOfDerivatives; l++) {
    			R3D[quadIndex][basisIndex][l].resize(numberOfDerivatives-l+1);
    			for (unsigned int m=0; m<=numberOfDerivatives-l; m++) {
      			R3D[quadIndex][basisIndex][l][m].resize(numberOfDerivatives-l-m+1, 0.0);
      			for (unsigned int n=0; n<=numberOfDerivatives-l-m; n++) {
							double v = 0.0; 
							for (unsigned int i=0; i<=l; i++) {
	  						double v1=0.0;
	  						for (unsigned int j=0; j<=m; j++) {
	    						double v2=0.0;
	    						for (unsigned int k=0; k<=n; k++) {
	      						if ((i+j+k)==0) 	v  = A[basisIndex][l][m][n];  //i=j=k=0
	      						else			v2+= Bin(n,k)*W[i][j][k]*R3D[quadIndex][basisIndex][l-i][m-j][n-k];
	    						}
	    						v1 += Bin(m,j)*v2;
	  						}
	  						v -= Bin(l,i)*v1;
							}
							R3D[quadIndex][basisIndex][l][m][n]=v/W[0][0][0];
     			 	}
    		 	}
  		 	}
		 	}
   	}
			
    //evaluate jacobian and compute quad point positions in physical space
    dealii::FullMatrix<double> Jacobian(dim, dim); std::vector<std::vector<double> >  basisdx;  
    dealii::Tensor<3,dim> dJ; std::vector<std::vector<std::vector<double> > > basisddx;
    basisdx.resize(controlpoints_per_cell);  
    if (numberOfDerivatives==2) basisddx.resize(controlpoints_per_cell);
    quadPointLocations[quadIndex].resize(dim, 0.0);
    //compute N, N,xi N,xixi
    for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
			//compute N and quadPointLocations
			double shapeValue;
			if (dim==1) 		shapeValue=R1D[quadIndex][basisIndex][0];
			else if (dim==2) 	shapeValue=R2D[quadIndex][basisIndex][0][0];
			else			 	shapeValue=R3D[quadIndex][basisIndex][0][0][0];
			basisValue[quadIndex][basisIndex]=shapeValue;
			for (unsigned int i=0; i<dim; i++){
  			quadPointLocations[quadIndex][i]+=shapeValue*cell->controlPoints[basisIndex]->coords[i];
			}
			//compute N,xi N,xixi
			basisdx[basisIndex].resize(dim, 0.0);
			if (numberOfDerivatives==2) basisddx[basisIndex].resize(dim);
			for (unsigned int I=0; I<dim; I++){
  			double dx;
  			if (dim==1) 		dx=R1D[quadIndex][basisIndex][1];
  			else if (dim==2) {
    			if (I==0)		dx=R2D[quadIndex][basisIndex][1][0];
    			else 	  		dx=R2D[quadIndex][basisIndex][0][1];
  			}
  			else if (dim==3) {
    			if (I==0) 		dx=R3D[quadIndex][basisIndex][1][0][0];
    			else if (I==1)  dx=R3D[quadIndex][basisIndex][0][1][0];
    			else			dx=R3D[quadIndex][basisIndex][0][0][1];
  			}
  			basisdx[basisIndex][I]=dx;
  			if (numberOfDerivatives==2){
    			basisddx[basisIndex][I].resize(dim, 0.0);
    			for (unsigned int J=0; J<dim; J++){
      			double ddx;
      			if (dim==1) 		ddx=R1D[quadIndex][basisIndex][2];
      			else if (dim==2) {
							if (I==0){
	  						if (J==0)	ddx=R2D[quadIndex][basisIndex][2][0];
	  						else		ddx=R2D[quadIndex][basisIndex][1][1];
							}
							else{
	  						if (J==0)	ddx=R2D[quadIndex][basisIndex][1][1];
	  						else		ddx=R2D[quadIndex][basisIndex][0][2];
							}
      			}
      			else if (dim==3) {
							if (I==0){
	  						if (J==0)		ddx=R3D[quadIndex][basisIndex][2][0][0];
	  						else if (J==1)	ddx=R3D[quadIndex][basisIndex][1][1][0];
	  						else			ddx=R3D[quadIndex][basisIndex][1][0][1];
							}
							else if (I==1){
	  						if (J==0)		ddx=R3D[quadIndex][basisIndex][1][1][0];
	  						else if (J==1)	ddx=R3D[quadIndex][basisIndex][0][2][0];
	  						else			ddx=R3D[quadIndex][basisIndex][0][1][1];
							}
							else{
	  						if (J==0)		ddx=R3D[quadIndex][basisIndex][1][0][1];
	  						else if (J==1)	ddx=R3D[quadIndex][basisIndex][0][1][1];
	  						else			ddx=R3D[quadIndex][basisIndex][0][0][2];											
							}
      			}
      			basisddx[basisIndex][I][J]=ddx;
    			}
  			}
			}
    }
			
    //evaluate J, dJ
    for (unsigned int i=0; i<dim; i++){
			for (unsigned int j=0; j<dim; j++){
				Jacobian[i][j]=0.0;
 			 	for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
   				Jacobian[i][j]+=basisdx[basisIndex][j]*cell->controlPoints[basisIndex]->coords[i];
  			}

  			if (numberOfDerivatives==2){
    			for (unsigned int k=0; k<dim; k++){
      			dJ[i][j][k]=0;
      			for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
							dJ[i][j][k]+=basisddx[basisIndex][j][k]*cell->controlPoints[basisIndex]->coords[i];
      			}
    			}
  			}
			}
    }
    dealii::FullMatrix<double> invJ(Jacobian); invJ.gauss_jordan(); ///inverse jacobian matrix, the "full size matrix"
    
    //find det(jacobian). But reduce the jacobian to (dim-1)x(dim-1) for quadrature points on the faces
    if (tempFaceQuadMap.count(quadIndex)>0){
			for (unsigned int i=0; i<dim; i++){
  			Jacobian[i][tempFaceQuadMap[quadIndex]]=Jacobian[tempFaceQuadMap[quadIndex]][i]=(i==tempFaceQuadMap[quadIndex]);
			}
    }
    //
    if (tempPlaneQuadMap.count(quadIndex)>0){
			for (unsigned int i=0; i<dim; i++){
  			Jacobian[i][tempPlaneQuadMap[quadIndex]]=Jacobian[tempPlaneQuadMap[quadIndex]][i]=(i==tempPlaneQuadMap[quadIndex]);
			}
    }
    //not line quad_point
    if(tempLineQuadMap.count(quadIndex)==0){
      jacobian[quadIndex]=std::abs(Jacobian.determinant())*(quadPoints[quadIndex][dim]); 
    //if (jacobian[quadIndex]<1.0e-15){printf("jacobian[quadIndex]<1.0e-15.\n"); exit(-1);}

      double parametricElementVolume=1.0;
      for (unsigned int i=0; i<dim; i++) {
  			if(tempFaceQuadMap.count(quadIndex)==0 and tempPlaneQuadMap.count(quadIndex)==0) {parametricElementVolume*=(cell->endKnots[i][1]-cell->endKnots[i][0])/2.0;}
  			if ( (tempFaceQuadMap.count(quadIndex)>0) and (i!=tempFaceQuadMap[quadIndex]) ) {parametricElementVolume*=(cell->endKnots[i][1]-cell->endKnots[i][0])/2.0;} 
        if ( (tempPlaneQuadMap.count(quadIndex)>0) and (i!=tempPlaneQuadMap[quadIndex]) ) {parametricElementVolume*=(cell->endKnots[i][1]-cell->endKnots[i][0])/2.0;}
       //in IGA here it's cell's level calculation, integral domain is (endKnots[i][1],endKnots[i][0]). Quadradual integration domian is from (-1,1) 
      }
      jacobian[quadIndex]*=parametricElementVolume;
    }
    if(tempLineQuadMap.count(quadIndex)>0){
			int lineDim=tempLineQuadMap[quadIndex];
			jacobian[quadIndex]=std::abs(Jacobian[lineDim][lineDim])*(quadPoints[quadIndex][dim])*(cell->endKnots[lineDim][1]-cell->endKnots[lineDim][0])/2.0;
    }
			
    //evaluate N,x from N,xi
    basisDx[quadIndex].resize(controlpoints_per_cell); 
    if (numberOfDerivatives==2) basisDDx[quadIndex].resize(controlpoints_per_cell);
    for (basisIndex=0;  basisIndex<controlpoints_per_cell; basisIndex++){
			basisDx[quadIndex][basisIndex].resize(dim, 0.0); 
			if (numberOfDerivatives==2) basisDDx[quadIndex][basisIndex].resize(dim);
			for (unsigned int I=0; I<dim; I++){
				double tempVal1=0.0;
				for (unsigned int i=0; i<dim; i++){
					tempVal1+=basisdx[basisIndex][i]*invJ[i][I];
				}
				basisDx[quadIndex][basisIndex][I]=tempVal1;
  			if (numberOfDerivatives==2) {
  				basisDDx[quadIndex][basisIndex][I].resize(dim, 0.0);
   			 for (unsigned int J=0; J<dim; J++){
   				double tempVal1=0.0 ,tempVal2=0.0;
    			for (unsigned int i=0; i<dim; i++){
						for (unsigned int j=0; j<dim; j++){
	  					tempVal1+= basisddx[basisIndex][i][j]*invJ[i][I]*invJ[j][J];
	  					for (unsigned int K=0; K<dim; K++){
	    					for (unsigned int l=0; l<dim; l++){
	     					 tempVal2+= basisdx[basisIndex][i]*invJ[i][K]*dJ[K][l][j]*invJ[l][I]*invJ[j][J];										
	   						}
	 						}
				 	 	}
     	 		}
     	 	 	basisDDx[quadIndex][basisIndex][I][J]=tempVal1-tempVal2;
     			//printf ("%12.6e %12.6e\n", tempVal1, tempVal2);
    			}
  			}
			}
		}
	}
}

/*
*following "return" functions are for users to access the values
*
*/
//coords of quadPoints
template <int dim>
std::vector<double>& IGAValues<dim>::quadraturePoints_coords( unsigned int q){
  return QuadPoints_coords[q];
}

//JxW
template <int dim>
double IGAValues<dim>::JxW(unsigned int q){
  return jacobian[q]; 
}

template <int dim>
double IGAValues<dim>::JxW_face(unsigned int q, unsigned int faceID){
  return jacobian[quadratureMap[faceID][q]]; 
}

template <int dim>
double IGAValues<dim>::JxW_line(unsigned int q, unsigned int lineID){
  return jacobian[LinequadratureMap[lineID][q]]; 
}

template <int dim>
double IGAValues<dim>::JxW_plane(unsigned int q, unsigned int planeID){
  return jacobian[PlanequadratureMap[planeID][q]]; 
}

template <int dim>
//shape_value
double IGAValues<dim>::shape_value(unsigned int i, unsigned int q){
  return basisValue[q][system_to_controlpoint_index(i)];
}

template <int dim>
double IGAValues<dim>::shape_value_face(unsigned int i, unsigned int q, unsigned int faceID){
  return basisValue[quadratureMap[faceID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
double IGAValues<dim>::shape_value_line(unsigned int i, unsigned int q, unsigned int lineID){
  return basisValue[LinequadratureMap[lineID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
double IGAValues<dim>::shape_value_plane(unsigned int i, unsigned int q, unsigned int planeID){
  return basisValue[PlanequadratureMap[planeID][q]][system_to_controlpoint_index(i)];
}

//shape_grad
template <int dim>
std::vector<double>& IGAValues<dim>::shape_grad(unsigned int i, unsigned int q){
  return basisDx[q][system_to_controlpoint_index(i)];
}

template <int dim>
std::vector<double>& IGAValues<dim>::shape_grad_face(unsigned int i, unsigned int q, unsigned int faceID){
  return basisDx[quadratureMap[faceID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
std::vector<double>& IGAValues<dim>::shape_grad_line(unsigned int i, unsigned int q, unsigned int lineID){
  return basisDx[LinequadratureMap[lineID][q]][system_to_controlpoint_index(i)];
}  

template <int dim>
std::vector<double>& IGAValues<dim>::shape_grad_plane(unsigned int i, unsigned int q, unsigned int planeID){
  return basisDx[PlanequadratureMap[planeID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
//shape_grad_grad
std::vector<std::vector<double> >& IGAValues<dim>::shape_grad_grad(unsigned int i, unsigned int q){
  return basisDDx[q][system_to_controlpoint_index(i)];
}

template <int dim>
std::vector<std::vector<double> >& IGAValues<dim>::shape_grad_grad_face(unsigned int i, unsigned int q, unsigned int faceID){
  return basisDDx[quadratureMap[faceID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
std::vector<std::vector<double> >& IGAValues<dim>::shape_grad_grad_line(unsigned int i, unsigned int q, unsigned int lineID){
  return basisDDx[LinequadratureMap[lineID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
  std::vector<std::vector<double> >& IGAValues<dim>::shape_grad_grad_plane(unsigned int i, unsigned int q, unsigned int planeID){
  return basisDDx[PlanequadratureMap[planeID][q]][system_to_controlpoint_index(i)];
}

template <int dim>
unsigned int IGAValues<dim>::system_to_component_index(unsigned int dof){
  if (dof>=dofs_per_cell) {printf("Requested dof greater than available regular dof. Are you trying to access special DOF like lagrange multipliers\n"); exit(-1);}
  std::div_t divresult= std::div ((int) dof, (int) dofPerControlPoint);
  return divresult.rem;
}

template <int dim>
unsigned int IGAValues<dim>::system_to_controlpoint_index(unsigned int dof){
  if (dof>=dofs_per_cell) {printf("Requested dof greater than available regular dof. Are you trying to access special DOF like lagrange multipliers\n"); exit(-1);}
  std::div_t divresult= std::div ((int) dof, (int) dofPerControlPoint);
  return divresult.quot;
}

template class IGAValues<1>;
template class IGAValues<2>;
template class IGAValues<3>;