//#include "Eigen/Sparse"

// Function declarations

void sfbc(VectorXd* sf,
		int left_lower,int left_increment,int left_upper,
		int bottom_lower,int bottom_increment,int bottom_upper,
		int right_lower,int right_increment,int right_upper,
		int top_lower,int top_increment,int top_upper);
		


void sjbc(SparseMatrix<double>* sj, int nx, int ny,
		int left_lower,int left_increment,int left_upper,
		int bottom_lower,int bottom_increment,int bottom_upper,
		int right_lower,int right_increment,int right_upper,
		int top_lower,int top_increment,int top_upper);
		
void zero_specific_row_sj(SparseMatrix<double>* sj, int row, int nx, int ny);
		
// Function structures

void sfbc(VectorXd* sf,
		int left_lower,int left_increment,int left_upper,
		int bottom_lower,int bottom_increment,int bottom_upper,
		int right_lower,int right_increment,int right_upper,
		int top_lower,int top_increment,int top_upper)
{
	for (unsigned i=left_lower;i<=left_upper;i=i+left_increment) {
		sf[0](i)=0.0;
	}
	for (unsigned i=bottom_lower;i<=bottom_upper;i=i+bottom_increment) {
		sf[0](i)=0.0;
	}
	for (unsigned i=top_lower;i<=top_upper;i=i+top_increment) {
		sf[0](i)=0.0;
	}
	for (unsigned i=right_lower;i<=right_upper;i=i+right_increment) {
		sf[0](i)=0.0;
	}
}

void sjbc(SparseMatrix<double>* sj, int nx, int ny,
		int left_lower,int left_increment,int left_upper,
		int bottom_lower,int bottom_increment,int bottom_upper,
		int right_lower,int right_increment,int right_upper,
		int top_lower,int top_increment,int top_upper) 
{
	for (unsigned i=left_lower;i<=left_upper;i=i+left_increment) {
		zero_specific_row_sj(sj,i,nx,ny);
		sj[0].coeffRef(i,i)=1.0;
	}
	for (unsigned i=bottom_lower;i<=bottom_upper;i=i+bottom_increment) {
		zero_specific_row_sj(sj,i,nx,ny);
		sj[0].coeffRef(i,i)=1.0;
	}
	for (unsigned i=top_lower;i<=top_upper;i=i+top_increment) {
		zero_specific_row_sj(sj,i,nx,ny);
		sj[0].coeffRef(i,i)=1.0;
	}
	for (unsigned i=right_lower;i<=right_upper;i=i+right_increment) {
		zero_specific_row_sj(sj,i,nx,ny);
		sj[0].coeffRef(i,i)=1.0;
	}
}

void zero_specific_row_sj(SparseMatrix<double>* sj, int row, int nx, int ny) {
	// CHANGE THE FORMATE TO MAKE IT LOOK NICER LATER
	const int remainder=row % 4;
	const int element=(row-remainder)/4;
	int left, middle, right;
	if (element<ny) { //Left
		if (element==0) { // Bottom left
			middle=0;
			right=4*ny;
			for (unsigned int i=0;i<8;i++) {
				sj[0].coeffRef(row,middle+i)=0.0;
				sj[0].coeffRef(row,right+i)=0.0;
			}
		}
		else if (element==ny-1) { // Top left
			middle=4*(element-1);
			right=middle+4*ny;
			for (unsigned int i=0;i<8;i++) {
				sj[0].coeffRef(row,middle+i)=0.0;
				sj[0].coeffRef(row,right+i)=0.0;
			}
		}
		else { // Middle left
			middle=4*(element-1);
			right=middle+4*ny;
			for (unsigned int i=0;i<12;i++) {
				sj[0].coeffRef(row,middle+i)=0.0;
				sj[0].coeffRef(row,right+i)=0.0;
			}
		}
	}
	else if (element % ny == 0) { //Bottom
		if (element == ny*(nx-1)) { //Bottom right
			left=4*ny*(nx-2);
			middle=4*ny*(nx-1);
			for (unsigned int i=0;i<8;i++) {
				sj[0].coeffRef(row,left+i)=0.0;
				sj[0].coeffRef(row,middle+i)=0.0;
			}
		}
		else { //Middle bottom
			middle=4*element;
			left=middle-4*ny;
			right=middle+4*ny;
			for (unsigned int i=0;i<8;i++) {
				sj[0].coeffRef(row,left+i)=0.0;
				sj[0].coeffRef(row,middle+i)=0.0;
				sj[0].coeffRef(row,right+i)=0.0;
			}
		}
	}
	else if ((element+1) % ny == 0) { // Top
		if (element==nx*ny-1) { //Top right
			middle=4*(element-1);
			left=middle-4*ny;
			for (unsigned int i=0;i<8;i++) {
				sj[0].coeffRef(row,left+i)=0.0;
				sj[0].coeffRef(row,middle+i)=0.0;
			}
		}
		else { // Middle top
			middle=4*(element-1);
			left=middle-4*ny;
			right=middle+4*ny;
			for (unsigned int i=0;i<8;i++) {
				sj[0].coeffRef(row,left+i)=0.0;
				sj[0].coeffRef(row,middle+i)=0.0;
				sj[0].coeffRef(row,right+i)=0.0;
			}
		}
	}
	else if (element > ny*(nx-1) ) { //Middle Right
		middle=4*(element-1);
		left=middle-4*ny;
		for (unsigned int i=0;i<12;i++) {
			sj[0].coeffRef(row,left+i)=0.0;
			sj[0].coeffRef(row,middle+i)=0.0;
		}
	}
	else { // Middle
		middle=4*(element-1);
		left=middle-4*ny;
		right=middle+4*ny;
		for (unsigned int i=0;i<12;i++) {
			sj[0].coeffRef(row,left+i)=0.0;
			sj[0].coeffRef(row,middle+i)=0.0;
			sj[0].coeffRef(row,right+i)=0.0;
		}
	}
}
