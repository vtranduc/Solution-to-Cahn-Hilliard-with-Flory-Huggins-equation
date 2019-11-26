#include <iostream>
using namespace std;

#include <cmath>
#include <cstdlib>
#include <fstream>
//#include <typeinfo> // For printing types. To be deleted!
//#include <iomanip>

#include "Eigen/Sparse"
#include "Eigen/CholmodSupport"
using namespace Eigen;

#include "basis.h"
#include "accessory.h"
#include "ic_bc.h"


//#include "Eigen/Sparse"


// Global variables are specified

const double diff = 6.0*pow(10.,2.0);
const int nex = 10, ney = 10; // 800x800 has resulted in Segmentation failt (core dumped)
// BE SURE TO CHANGE LINE 521 (APPROXIMATELY), SINCE NUMBER OF DIAGONALS WILL VARY!

const int nx = nex + 1, ny = ney + 1;

//List internal functions

// Functions used for applying BC=================================================================

const int left_lower=1;
const int left_increment=4;
const int left_upper=1+4*ney;
const int bottom_lower=2;
const int bottom_increment=4*ny;
const int bottom_upper=2+4*nex*ny;
const int top_lower=2+4*ney;
const int top_increment=4*ny;
const int top_upper=2+4*ney+4*nex*ny;
const int right_lower=1+4*nex*ny;
const int right_increment=4;
const int right_upper=1+4*ney+4*nex*ny;

void cIC(double *c) {
	for (unsigned i=1;i<=left_upper;i=i+4) {
		c[i]=0.0;
	}
	for (unsigned i=2;i<=bottom_upper;i=i+bottom_increment) {
		c[i]=0.0;
	}
	for (unsigned i=top_lower;i<=top_upper;i=i+top_increment) {
		c[i]=0.0;
	}
	for (unsigned i=right_lower;i<=right_upper;i=i+4) {
		c[i]=0.0;
	}
}

void sfBC(VectorXd *sf) {
	for (unsigned i=1;i<=left_upper;i=i+4) {
		sf[0](i)=0.0;
	}
	for (unsigned i=2;i<=bottom_upper;i=i+bottom_increment) {
		sf[0](i)=0.0;
	}
	for (unsigned i=top_lower;i<=top_upper;i=i+top_increment) {
		sf[0](i)=0.0;
	}
	for (unsigned i=right_lower;i<=right_upper;i=i+4) {
		sf[0](i)=0.0;
	}
}

void sjBC(SparseMatrix<double> *sj, int rows, int cols) {
	for (unsigned i=1;i<=left_upper;i=i+4) {
		for (unsigned j=0; j<cols; j++) {
			if (sj[0].coeffRef(i,j)!='\0' && sj[0].coeffRef(i,j)!=0.0) { // Perhaps unnessessary to check for NULL!!!
				sj[0].coeffRef(i,j)=0.0;
			}
		}
		sj[0].coeffRef(i,i)=1.0;
	}
	for (unsigned i=2;i<=bottom_upper;i=i+bottom_increment) {
		for (unsigned j=0; j<cols; j++) {
			if (sj[0].coeffRef(i,j)!='\0' && sj[0].coeffRef(i,j)!=0.0) {
				sj[0].coeffRef(i,j)=0.0;
			}
		}
		sj[0].coeffRef(i,i)=1.0;
	}
	for (unsigned i=top_lower;i<=top_upper;i=i+top_increment) {
		for (unsigned j=0; j<cols; j++) {
			if (sj[0].coeffRef(i,j)!='\0' && sj[0].coeffRef(i,j)!=0.0) {
				sj[0].coeffRef(i,j)=0.0;
			}
		}
		sj[0].coeffRef(i,i)=1.0;
	}
	for (unsigned i=right_lower;i<=right_upper;i=i+4) {
		for (unsigned j=0; j<cols; j++) {
			if (sj[0].coeffRef(i,j)!='\0' && sj[0].coeffRef(i,j)!=0.0) {
				sj[0].coeffRef(i,j)=0.0;
			}
		}
		sj[0].coeffRef(i,i)=1.0;
	}
}

// Other accessory functions ========================================================================



// ================================================================================================
// ================================================================================================
// ================================================================================================
// =============================== Start the main function ========================================
// ================================================================================================
// ================================================================================================
// ================================================================================================

int main() {

	// Start counting time =======================================
	
	clock_t t0=clock();

	// Do draft testing here======================================
	
	
	cout << "End of test ================================" << endl;
	
	//===========================================================
	
	cout << "Hello World!" << endl;
	
	const double w[3]={0.27778, 0.4444, 0.27778};
	const double gp[3]={0.1127, 0.5, 0.8873};
	
	const int ne = nex * ney;
	const int n = nx * ny;
	
	const int nfour = 4 * n;
	double time = 0.0;
	
	int k, l;
	
	//------------------------------------------------
	
	double dto = 1.0 * pow(10.,-5.), dt = 1.0 * pow(10.,-5.);
	double max_time = 1.2 * pow(10.,-4.);
	
	double dt_catchup = 5.0 * pow(10.,-6.);
	double catchup_pace = 1.5;
	
	//--------------------------------------------------
	
	// Positions of grid point========================================================================
	
	double x[n];
	double y[n];
	
	for (unsigned int i = 0; i < n; i++) {
		x[i] = (1.0/nex) * floor(i/ny);
	}
	
	double z1;
	for (unsigned int i = 0; i < ny; i++) {
		z1 = (1.0/ney)*i;
		for (unsigned int j = 0; j < nx; j++) {
			k = i+j*ny;
			y[k] = z1;
		}
	}
	
	cout << "====================================================" << endl;
	
	// Create matrices that help to characterize each element uniquely=================================
	
	int nop[ne][16];
	int nopm[ne][4];
	
	for (unsigned int i=0; i<ne; i++) {
		nop[i][0]=4*(i/ney)+4*i+1;
		nop[i][8]=nop[i][0]+4*ny;
		for (unsigned int j=0; j<7; j++) {
			k=1+j;
			l=9+j;
			nop[i][k]=nop[i][0]+k;
			nop[i][l]=nop[i][8]+k;
		}
	}
	
	for (unsigned int i=0; i<ne; i++) {
		nopm[i][0]=i/ney+i+1;
		nopm[i][1]=nopm[i][0]+1;
		nopm[i][2]=nopm[i][0]+ny;
		nopm[i][3]=nopm[i][2]+1;
	}
	
	// Randomize initial conditions ========================================================================
	
	double co[nfour];
	double LO=0.7-0.01, HI=0.7+0.01;
	for (unsigned int i=0; i<nfour; i++) {
		co[i] = LO + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(HI-LO)));
	}

	// Take the first time step =============================================================================
	
	time=time+dto;
	double c[nfour];	
	equateArrays(c, co, nfour);
	
	// Apply boundary conditions ===========================================
	
	cIC(&c[0]);
	
	// Export the data of initial conds here===================================================================
	
	// Prepare main loop ====================================================================================
	
	double cp[nfour];
	double coo[nfour];
	
	double chi=1.1666666666667;
	
	double phi[16], phix[16], phiy[16], phixx[16], phiyy[16], phixy[16];
	double con, cono, conx, cony, conxx, conyy;
	double cont;
	
	int inop[16], inopm[4];
	
	double dx, dy;
	long double error;
	
	//cout.precision(17);
	
	//-----------------------------------------------------------------------------

	/*
	double **sj; THIS IS STILL WORTH NOTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	sj = (double**) malloc(sizeof(double)*nfour);
	for (unsigned int i=0; i<nfour; i++) {
		sj[i] = (double*) malloc(sizeof(double)*nfour);
	}
	*/
	
	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	
	VectorXd c_(nfour), sf(nfour);
	SparseMatrix<double> sj(nfour,nfour);
	
	cout << "Initializing" << endl;
	sj.reserve(VectorXi::Constant(nfour,36));
	testIni(&sj,nx, ny);
	
	//sj.insert(5,7)=0.0;
	
	//cout << "Compressing " << sj.coeffRef(nfour-1,0) << endl;
	sj.makeCompressed();
	cout << "Current nonZero: " << sj.nonZeros() << endl;
	//cout << sj << endl;
	cout << "Generate handler" << endl;
	BiCGSTAB<SparseMatrix<double> > handler;
	
	//CholmodBase<SparseMatrix<double> > handler;
	
	//SparseLU<SparseMatrix<double> > handler;
	
	//IncompleteLUT<SparseMatrix<double>,IncompleteLUT<double> > handler;
	
	//BiCGSTAB< SparseMatrix<double>,IncompleteLUT<double> > handler;
	
	//handler.setMaxIterations(3);
	
	handler.setTolerance(0.000001);
	
	handler.analyzePattern(sj);
	
	cout << "Enter main loop" << endl;
	//cout << sj << endl;
	
	//------------------------------------------------------------------------------
	
	int jk;
	
	// Enter Newton-Raphson & Time predictor loop ===============================================
	
	for (int rahmat=0;rahmat<1;rahmat++) { //REPLACE IT WITH FOREVER LOOP LATER
	
	// Prediction & stores previous results in case of error ====================================================
	
		for (unsigned int i=0; i<nfour; i++) {
			cp[i]=c[i]+dt*((c[i]-co[i])/dto);
		}
		equateArrays(coo,co,nfour);
		equateArrays(co,c,nfour);
		equateArrays(c,cp,nfour);
		
	// Apply BC to newly predicted c =============================================================================
		
		cIC(c); // THIS IS PERHAPS NOT NECESSARY IF WE APPLY BC TO COO BEFORE MAIN LOOP STARTS!!!!!
		
	// Prepare to enter Newton-Raphson loop =========================================================================
	
		jk=0;
	
	// Enter Newton-Raphson loop ==================================================================================
	
		for (int saeedi=0; saeedi<2; saeedi++) { // REPLACE IT WITH FOREVER LOOP LATER
		
	// Reset the residual vector and Jacobian matrix ===============================================================
		
			//resetArray(&sf[0], nfour);
			resetsf(&sf, nfour);
			resetsj(&sj,nx,ny); //THIS LINE TAKES ABOUT 4 SECONDS!!! PERHAPS IT'S BETTER TO USE MALLOC!!!

			//cout << sj<< endl;
			
	// Start looping each element ==========================================================================
	cout << " start filling " << endl;
			for (unsigned int i=0; i<ne; i++) {
				//cout << "Working on " << i << "th element" << endl;
	// Characterize the element ============================================================================
	
				for (unsigned int j=0; j<4; j++) {
					inopm[j]=nopm[i][j]-1; //inopm is used solely for index, so we decrement by one
				}
				for (unsigned int j=0; j<16; j++) {
					inop[j]=nop[i][j]-1; // we decrement by one for the same reason as above
				}
				dx = x[inopm[2]] - x[inopm[0]];
				dy = y[inopm[1]] - y[inopm[0]];
				
	// Shift the weights of two basis functions ===========================================================
	
				for (unsigned int j=0; j<=2; j++) {
					for (unsigned int k=0; k<=2; k++) {
						
	// Characterize the current element's basis functions according to the weight =========================

						fPHI(phi,gp[j],gp[k],dx,dy);
						fPHIX(phix,gp[j],gp[k],dx,dy);
						fPHIY(phiy,gp[j],gp[k],dx,dy);
						fPHIXX(phixx,gp[j],gp[k],dx,dy);
						fPHIYY(phiyy,gp[j],gp[k],dx,dy);
						fPHIXY(phixy,gp[j],gp[k],dx,dy);
						
	// Determine previous and current absolute values, and associated slopes at current particular point ==
	
						con=0.0;
						cono=0.0;
						conx=0.0;
						cony=0.0;
						conxx=0.0;
						conyy=0.0;
						
						for (unsigned int intr=0; intr<16; intr++) {
							con = con + c[inop[intr]] * phi[intr];
		                    cono = cono + co[inop[intr]] * phi[intr];
		                    conx = conx + c[inop[intr]] * phix[intr];
		                    cony = cony + c[inop[intr]] * phiy[intr];
		                    conxx = conxx + c[inop[intr]] * phixx[intr];
		                    conyy = conyy + c[inop[intr]] * phiyy[intr];
						}
						
	// Change in concentration over time step is estimated ================================================
	
						cont = (con - cono) / dto; // In original code, dt is used instead
						
	// Fill in residual vector ============================================================================
	
						for (unsigned int l=0; l<16; l++) {
							sf(inop[l])=sf(inop[l])-w[j]*w[k]*dx*dy*
								(cont*phi[l]-diff*phi[l]*(-1.0/(pow(con,2.0))
								+(1./10.)*pow((1-con),-2.0))*(pow(conx,2.0)+
								pow(cony,2.0))-diff*phi[l]*
								(1/con+(1./10.)*(1./(1.-con))-2*chi)*
								(conxx+conyy)+(conxx+conyy)*
								(phixx[l]+phiyy[l]));
							
								
	// Fill in Jacobian matrix ============================================================================
							
							for (unsigned int m=0; m<16; m++) {
								sj.coeffRef(inop[l],inop[m])=sj.coeffRef(inop[l],inop[m])+w[j]
		                            *w[k]*dx*dy*
		                            (phi[l]*phi[m]/dt-diff*phi[l]*((2*phi[m]
		                            /(pow(con,3.0))+2*phi[m]/(10*pow((1-con),3)))*
		                            (pow(conx,2.0)+pow(cony,2.0))+(-1/pow(con,2.0)+1/(10*pow((1-con),2.0)))
		                            *(2*conx*phi[m]+2*cony*phi[m]))-diff*phi[l]
		                            *((-phi[m]/pow(con,2.0)+phi[m]/(10*pow((1-con),2)))*
		                            (conxx+conyy)+(1/con+1/(10*(1-con))-2*chi)
		                            *(phixx[m]+phiyy[m]))+(phixx[m]+phiyy[m])
									*(phixx[l]+phiyy[l]));
									
							}
						}
					}
				}
			}
			
	// Apply boundary conditions ===========================================================================
			cout << "Current nonZero: " << sj.nonZeros() << endl;
			cout << "Applying BC..." << endl;
			sfbc(&sf,
				left_lower,left_increment,left_upper,
				bottom_lower,bottom_increment,bottom_upper,
				right_lower,right_increment,right_upper,
				top_lower,top_increment,top_upper);
			
			sjbc(&sj,nx, ny,
				left_lower,left_increment,left_upper,
				bottom_lower,bottom_increment,bottom_upper,
				right_lower,right_increment,right_upper,
				top_lower,top_increment,top_upper);
				
			//sjBC(&sj, nfour, nfour);
			cout << "Current nonZero: " << sj.nonZeros() << endl;
			//cout << sf<< endl;
			//cout << sj << endl;
			
	// Solve linear system =================================================================================
	///*
			cout << "Solving matrix" << endl;
			cout << "Factorize" << endl;
			handler.factorize(sj);
			cout << "solve" << endl;
			c_=handler.solve(sf);
			
			cout << "#iterations:     " << handler.iterations() << endl;
			cout << "estimated error: " << handler.error()      << endl;
			
	// Update the solution =================================================================================
	
			//cout << sj << endl;
			cout << "Update the solution" << endl;
			for (unsigned int i=0; i<nfour; i++) {
				c[i]=c[i]+c_(i);
			}
			
	// Evaluate the error ==================================================================================
	
			//error=evalError(&c_,nfour);
			//cout << setprecision(13) << fixed  <<"Error for this loop is: " << error << endl;
			error = 0.0;
			for (unsigned int i=0; i<nfour; i++) {
				error=error+pow(c_(i),2.0);
			}
			error=sqrt(error);
			cout << "Error: " << error << endl;
			//*/
		}
	}
	
	// Conclude the program with elapsed time ==============================================================

	clock_t t1=clock();
	cout << "Elapsed time: " << double(t1-t0)/CLOCKS_PER_SEC << " seconds" << endl;
	
	return 0;
}
