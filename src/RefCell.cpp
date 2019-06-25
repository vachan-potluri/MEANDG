#include "RefCell.h"


//_________________________________________________Base class__________________________________________________//
// Reference cell

//constructor
RefCell :: RefCell(){
};

//destructor
RefCell :: ~RefCell(){
};

int RefCell :: getNoOfPoints()const{
	return this->noOfPoints;
};

void RefCell :: setNoOfPoints(int pts){
	this->noOfPoints = pts;
};

int RefCell :: getNoOfFaces()const{
	return this->noOfFaces;
};

void RefCell :: setNoOfFaces(int fcs){
	this->noOfFaces = fcs;
};

int RefCell :: getNoOfDOFPoints()const{
	return this->nDOF;
};

void RefCell :: setNoOfDOFPoints(int dofs){
	this->nDOF = dofs;
};

void RefCell::setIntegrationType(IntFlag::intflag intFlag){
	assert(intFlag == IntFlag::inexact or intFlag == IntFlag::exact && "Set integration type to inexact (intFlag 0) or exact (intFlag 1)");
	this->intFlag = intFlag;
};


TensorO1<double>* RefCell::getQuadWeights(){
	return &this->w_q;
};


void RefCell::init(int order, IntFlag::intflag intFlag, CellType::cellType celltype) {
	// order: polynomial order. 
	// intFlag: integration type. 0:inexact, 1:exact
	assert(intFlag == 0 or intFlag == 1 && "Set integration type to inexact (intFlag 0) or exact (intFlag 1)");

	// Set the cell type
	this->celltype = celltype; 
	// Set the order
	this->order = order;
	// Set the integration flag (exact vs inexact)
	this->intFlag = intFlag;

	assert (celltype == CellType::Hex or celltype == CellType::Tet or celltype == CellType::Prism or celltype == CellType::Pyramid);

	if (celltype == CellType::Hex){
		this->noOfPoints = 8;
		this->noOfFaces = 6;
		this->nDOF = pow((this->order + 1),3); //coded for a 3d goemetry only.

		int K = pow((this->order + 1 + this->intFlag),3);

		//setting up the cell Jacobian at quadrature points. For a refCell, all det(J) = 1.0
		J.setSize(K);

		for (int i=0; i<K; i++){
			J.setValue(i, 1.0);
		};
	}

	else if (celltype == CellType::Tet){
	}
	else if (celltype == CellType::Prism){
	}
	else if (celltype == CellType::Pyramid){
	};

	this->generateQuadratureWeights();
	this->generateLagMatrix();
	this->generateGradLagMatrix();
}


CellType::cellType RefCell::getCellType(){
	return this->celltype;
};

void RefCell::generateQuadratureWeights(){
	// Fills array w_q
	FunctionalSpace F(this->order, this->intFlag);
	if (celltype == CellType::Hex){
		int K = this->order + 1 + this->intFlag;
		int nQuad = K*K*K; //total number of quadrature points.
		this->nQuad = nQuad;

		w_q.setSize(nQuad);
		TensorO1<double> x_dummy(nQuad); 
		TensorO1<double> y_dummy(nQuad); 
		TensorO1<double> z_dummy(nQuad);

		F.LGLRootsAndWeights3D(K, K, K, &x_dummy, &y_dummy, &z_dummy, &this->w_q);
			
	}
	else if (celltype == CellType::Tet){
	}
	else if (celltype == CellType::Prism){
	}
	else if (celltype == CellType::Pyramid){
	}
};



void RefCell::generateLagMatrix(){
	// LagMatrix is a matrix of the Lagrange polynomials computed at the quadrature points.
	// It has dim: nQuad x nDOF

	FunctionalSpace F(this->order, this->intFlag);
	if (celltype == CellType::Hex){
		F.generateLagMatrix3DTensor(&this->LagMatrix);
	}
	else if (celltype == CellType::Tet){
		//insert function for generating Lagrange matrix for tet cells
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Prism){
		//insert function for generating Lagrange matrix for prism
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Pyramid){
		//insert function for generating Lagrange matrix for pyramid
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
};



void RefCell::generateGradLagMatrix(){
	// LagDrstMatrix is a matrix of the gradient of Lagrange polynomials in r, s, t direction computed at the quadrature points.
	// It has dim: nQuad x nDOF 

	FunctionalSpace F(this->order, this->intFlag);
	if (celltype == CellType::Hex){
		F.generateGradLagMatrix3DTensor(&this->LagDrMatrix, &this->LagDsMatrix, &this->LagDtMatrix);
	}
	else if (celltype == CellType::Tet){
		//insert function for generating Lagrange matrix for tet
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Prism){
		//insert function for generating Lagrange matrix for prism
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
	else if (celltype == CellType::Pyramid){
		//insert function for generating Lagrange matrix for pyramid
		//write correpsonding LagMatrix matrix generator function in functionalSpaces.cpp
	}
};


void RefCell::generateProjMatrix()
{
	// Set size
	int n_subcells = this->order+1;
	this->ProjMatrix.setSize(n_subcells*n_subcells*n_subcells,(this->order+1)*(this->order+1)*(this->order+1));
	this->ProjMatrix.setAll(0.0);

	// Projection matrix requires order/2+1 LGL quadrature
	// Get these points and weights
	TensorO1<double> q_points_x, q_points_y, q_points_z, q_weights; // q for quadrature
	int temp = this->order/2 + 1;
	FunctionalSpace Fn(temp,IntFlag::exact);
	Fn.LGLRootsAndWeights3D(temp,temp,temp,&q_points_x,&q_points_y,&q_points_z,&q_weights);

	// Get interpolation points for evaluating lagrange polynomials
	TensorO1<double> i_points, dummy;
	Fn.setOrder(this->order);
	Fn.LGLRootsAndWeights1D(&i_points,&dummy);

	int xid, yid, zid, j; // IDs for subcell indexing
	int i; // index for lagrange basis indexing
	int k; // index for quadrature
	double x, y, z; // current quadrature location mapped to refcell
	double curr_contrib; // current quadrature point contribution to P_ji
	TensorO1<double> f; // for evaluating lagrange polynomial
	f.setAll(0.0);
	for(xid = 0; xid<n_subcells; xid++){
		for(yid = 0; yid<n_subcells; yid++){
			for(zid=0; zid<n_subcells; zid++){
				j = zid + n_subcells*yid + n_subcells*n_subcells*xid;
				for(i=0; i<pow(Fn.getOrder()+1,3); i++){
					f.setValue(i,1.0); // set current polynomial coefficient to 1
					for(k=0; k<q_weights.getSize(); k++){
						x = 2.0*xid/n_subcells - 1.0 + (q_points_x.getValue(k)+1.0)/n_subcells;
						y = 2.0*yid/n_subcells - 1.0 + (q_points_y.getValue(k)+1.0)/n_subcells;
						z = 2.0*zid/n_subcells - 1.0 + (q_points_z.getValue(k)+1.0)/n_subcells;
						curr_contrib = q_weights.getValue(k)*Test::lagrangeInterpolation3D(x,y,z,&i_points,&f);
						this->ProjMatrix.addValue(j,i,curr_contrib);
					}
					f.setValue(i,0.0); // reset f to all-zeros
				}
			}
		}
	}
}


// Get pointer to the Lagrange matrix
TensorO2<double>* RefCell::getLagMatrix(){
	return &this->LagMatrix;
};


// Get pointer to the Lagrange derivative matrix in x direction (r in reference frame)
TensorO2<double>* RefCell::getLagDrMatrix(){
	return &this->LagDrMatrix;
};

// Get pointer to the Lagrange derivative matrix in y direction (s in reference frame)
TensorO2<double>* RefCell::getLagDsMatrix(){
	return &this->LagDsMatrix;
};

// Get pointer to the Lagrange derivative matrix in z direction (t in reference frame)
TensorO2<double>* RefCell::getLagDtMatrix(){
	return &this->LagDtMatrix;
};



// Other printing functions

void RefCell::printLagMatrix(){
	cout <<"\nLagMat: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
};


void RefCell::printGradLagMatrix(){
	cout <<"\nGradLagMat r: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagDrMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
	cout <<"\nGradLagMat s: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagDsMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
	cout <<"\nGradLagMat t: \n";
	for (int i=0; i<nQuad; i++){
		for (int j=0; j<nDOF; j++){
			cout << LagDtMatrix.getValue(i,j) << "\t";
		};
		cout << endl;
	};
};

void RefCell::print(){
	string celltypes[4] = {"Hex", "Tet", "Prism", "Pyramid"}; 

	cout << "\nThe cell type is: " << celltypes[this->celltype] << endl;
	cout << "\nTotal number of defining points are: " << this->noOfPoints << endl;
	cout << "\nTotal number of defining faces are: " << this->noOfFaces << endl;
	cout << "\nOrder of polynomial reconstruction: " << this->order << endl;
	cout << endl;
	cout << "Quadrature weights:  [ ";
	for (int i=0; i<nQuad; i++){
		cout << w_q.getValue(i) << "  " ;
	};
	cout << " ]\n";
};



//*******************************************************************************************************************//
// Tests:

namespace Test{
	void TestReferenceCells(){

		cout << "\nRunning tests on the Reference Cells.\n";

		for (Index i=1; i<3; i++){
			RefCell hexCell;
			hexCell.init(i, IntFlag::inexact, CellType::Hex);
			hexCell.generateQuadratureWeights();

			assert(hexCell.getNoOfPoints() == 8);
			assert(hexCell.getNoOfFaces() == 6);
			assert(hexCell.getNoOfDOFPoints() == pow((i+1),3));
			assert(Math::approx(hexCell.getQuadWeights()->getL1Norm(), 8.0));

		};

		cout << "Test::TestReferenceCells() passed.\n";
	};
};

