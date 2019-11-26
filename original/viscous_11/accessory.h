//Accessory functions

//#include <iostream>
//using namespace std;

//#include "Eigen/Sparse"

//using namespace Eigen;

void equateArrays(double *array1, double *array2, int length) {
	for (unsigned int i=0; i<length; i++) {
		array1[i]=array2[i];
	}
}

void resetsf(VectorXd* sf, int n) {
	for (unsigned int i=0;i<n;i++) {
		sf[0](i)=0.0;
	}
}

void initializesj(SparseMatrix<double>* sj, int rows, int cols) {
	for (unsigned int i=0; i<rows; i++) {
		for (unsigned int j=0; j<cols; j++) {
			sj[0].insert(i,j)=0.0;
		}
	}
}

void testIni(SparseMatrix<double>* sj, int nx, int ny) {
	// CHANGE THE FORMATE TO MAKE IT LOOK NICER LATER
	const int n=nx*ny;
	int _bound, bound_;
	int lower,middle, higher;
	for (unsigned int i=0; i<n; i++) {
		//cout << "modding " << i % ny << endl;
		if (i<ny) { //Left
			if (i==0) { //Bottom left corner
				higher=4*ny;
				for (unsigned int j=0; j<4;j++) {
					for (unsigned int k=0;k<8;k++) {
						sj[0].insert(j,k)=0.0;
						sj[0].insert(j,higher+k)=0.0;
					}
				}
			}
			else if (i==ny-1) { //Top left corner
				_bound=4*(ny-1);
				bound_=_bound+4;
				middle=_bound-4;
				higher=middle+4*ny;
				for (unsigned int j=_bound; j<bound_;j++) {
					for (unsigned int k=0;k<8;k++) {
						sj[0].insert(j,middle+k)=0.0;
						sj[0].insert(j,higher+k)=0.0;
					}
				}
			}
			else { //Middle left
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					//cout << j << " is the element ml" << n<< endl;
					lower=_bound-4;
					higher=bound_+4;
					for (unsigned int k=lower;k<higher; k++) {
						sj[0].insert(j,k)=0.0;
					}
					lower=_bound+4*(ny-1);
					higher=lower+12;
					for (unsigned int k=lower;k<higher; k++) {
						sj[0].insert(j,k)=0.0;
					}
				}
			}
		}
		
		else if (i % ny == 0) { //Bottom
			if (i == ny*(nx-1)) { //Bottom right
				_bound=4*i;
				bound_=_bound+4;
				//cout << "CHECK POINT" << endl;

				for (unsigned int j=_bound; j<bound_;j++) {
					lower=_bound-4*ny;
					higher=lower+8;
					for (unsigned int k=lower;k<higher; k++) {
						//cout << "up " << j << " " << k << " " << higher << endl;
						sj[0].insert(j,k)=0.0;
					}
					lower=_bound;
					higher=lower+8;
					for (unsigned int k=lower;k<higher; k++) {
						//cout << "down " << j << " " << k << endl;
						sj[0].insert(j,k)=0.0;
					}
				}
			}
			else { //Middle bottom
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					for (unsigned int k=0; k< 3;k++) {
						lower=_bound-4*ny+4*ny*k;
						higher=lower+8;
						for (unsigned int l=lower;l<higher;l++) {
							//cout << "up " << j << " " << l << " " << higher << endl;
							sj[0].insert(j,l)=0.0;
						}
					}
				}
			}
		}
		
		else if ((i+1) % ny == 0) { //Top
			if (i==nx*ny-1) { //Top right
				//cout << "CHECK POgsdfgsdfgsdINT" << endl;
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					lower=_bound-4*(ny+1);
					higher=lower+8;
					for (unsigned int k=lower; k<higher; k++) {
						sj[0].insert(j,k)=0.0;
					}
					lower=lower+4*ny;
					higher=lower+8;
					for (unsigned int k=lower; k<higher; k++) {
						sj[0].insert(j,k)=0.0;
					}
				}
			}
			else { // Middle top
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					for (unsigned int k=0; k< 3;k++) {
						lower=_bound-4*(ny+1)+4*ny*k;
						higher=lower+8;
						for (unsigned int l=lower;l<higher;l++) {
							sj[0].insert(j,l)=0.0;
						}
					}
				}
			}
		}
		
		else if (i > ny*(nx-1) ) { //Right
			_bound=4*i;
			bound_=_bound+4;
			for (unsigned int j=_bound; j<bound_;j++) {
				for (unsigned int k=0; k< 2;k++) {
					lower=_bound-4*(ny+1)+4*ny*k;
					higher=lower+12;
					for (unsigned int l=lower;l<higher;l++) {
						sj[0].insert(j,l)=0.0;
					}
				}
			}
		}
		else { // Middle
			_bound=4*i;
			bound_=_bound+4;
			for (unsigned int j=_bound; j<bound_;j++) {
				for (unsigned int k=0; k<3;k++) {
					lower=_bound-4*(ny+1)+4*ny*k;
					higher=lower+12;
					for (unsigned int l=lower;l<higher;l++) {
						sj[0].insert(j,l)=0.0;
					}
				}
			}
		}
	}
}

void resetsj(SparseMatrix<double>* sj, int nx, int ny) {
	// CHANGE THE FORMATE TO MAKE IT LOOK NICER LATER
	const int n=nx*ny;
	int _bound, bound_;
	int lower,middle, higher;
	for (unsigned int i=0; i<n; i++) {
		//cout << "modding " << i % ny << endl;
		if (i<ny) { //Left
			if (i==0) { //Bottom left corner
				higher=4*ny;
				for (unsigned int j=0; j<4;j++) {
					for (unsigned int k=0;k<8;k++) {
						sj[0].coeffRef(j,k)=0.0;
						sj[0].coeffRef(j,higher+k)=0.0;
					}
				}
			}
			else if (i==ny-1) { //Top left corner
				_bound=4*(ny-1);
				bound_=_bound+4;
				middle=_bound-4;
				higher=middle+4*ny;
				for (unsigned int j=_bound; j<bound_;j++) {
					for (unsigned int k=0;k<8;k++) {
						sj[0].coeffRef(j,middle+k)=0.0;
						sj[0].coeffRef(j,higher+k)=0.0;
					}
				}
			}
			else { //Middle left
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					//cout << j << " is the element ml" << n<< endl;
					lower=_bound-4;
					higher=bound_+4;
					for (unsigned int k=lower;k<higher; k++) {
						sj[0].coeffRef(j,k)=0.0;
					}
					lower=_bound+4*(ny-1);
					higher=lower+12;
					for (unsigned int k=lower;k<higher; k++) {
						sj[0].coeffRef(j,k)=0.0;
					}
				}
			}
		}
		
		else if (i % ny == 0) { //Bottom
			if (i == ny*(nx-1)) { //Bottom right
				_bound=4*i;
				bound_=_bound+4;
				//cout << "CHECK POINT" << endl;

				for (unsigned int j=_bound; j<bound_;j++) {
					lower=_bound-4*ny;
					higher=lower+8;
					for (unsigned int k=lower;k<higher; k++) {
						//cout << "up " << j << " " << k << " " << higher << endl;
						sj[0].coeffRef(j,k)=0.0;
					}
					lower=_bound;
					higher=lower+8;
					for (unsigned int k=lower;k<higher; k++) {
						//cout << "down " << j << " " << k << endl;
						sj[0].coeffRef(j,k)=0.0;
					}
				}
			}
			else { //Middle bottom
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					for (unsigned int k=0; k< 3;k++) {
						lower=_bound-4*ny+4*ny*k;
						higher=lower+8;
						for (unsigned int l=lower;l<higher;l++) {
							//cout << "up " << j << " " << l << " " << higher << endl;
							sj[0].coeffRef(j,l)=0.0;
						}
					}
				}
			}
		}
		
		else if ((i+1) % ny == 0) { //Top
			if (i==nx*ny-1) { //Top right
				//cout << "CHECK POgsdfgsdfgsdINT" << endl;
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					lower=_bound-4*(ny+1);
					higher=lower+8;
					for (unsigned int k=lower; k<higher; k++) {
						sj[0].coeffRef(j,k)=0.0;
					}
					lower=lower+4*ny;
					higher=lower+8;
					for (unsigned int k=lower; k<higher; k++) {
						sj[0].coeffRef(j,k)=0.0;
					}
				}
			}
			else { // Middle top
				_bound=4*i;
				bound_=_bound+4;
				for (unsigned int j=_bound; j<bound_;j++) {
					for (unsigned int k=0; k< 3;k++) {
						lower=_bound-4*(ny+1)+4*ny*k;
						higher=lower+8;
						for (unsigned int l=lower;l<higher;l++) {
							sj[0].coeffRef(j,l)=0.0;
						}
					}
				}
			}
		}
		
		else if (i > ny*(nx-1) ) { //Right
			_bound=4*i;
			bound_=_bound+4;
			for (unsigned int j=_bound; j<bound_;j++) {
				for (unsigned int k=0; k< 2;k++) {
					lower=_bound-4*(ny+1)+4*ny*k;
					higher=lower+12;
					for (unsigned int l=lower;l<higher;l++) {
						sj[0].coeffRef(j,l)=0.0;
					}
				}
			}
		}
		else { // Middle
			_bound=4*i;
			bound_=_bound+4;
			for (unsigned int j=_bound; j<bound_;j++) {
				for (unsigned int k=0; k<3;k++) {
					lower=_bound-4*(ny+1)+4*ny*k;
					higher=lower+12;
					for (unsigned int l=lower;l<higher;l++) {
						sj[0].coeffRef(j,l)=0.0;
					}
				}
			}
		}
	}
}
