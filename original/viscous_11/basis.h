//Defining basis function meshes

double f1(double z) {
return 1.-3.*pow(z,2.0)+2.*pow(z,3.0);
}
double f2(double z,double dz) {
return (z-2.*pow(z,2.0)+pow(z,3.0))*dz;
}
double f3(double z) {
return 3.*pow(z,2.0)-2.*pow(z,3.0);
}
double f4(double z,double dz) {
return (pow(z,3.0)-pow(z,2.0))*dz;
}
double f1z(double z,double dz) {
return (-6.*z+6.*pow(z,2.0))/dz;
}
double f2z(double z) {
return 1.-4.*z+3.*pow(z,2.0);
}
double f3z(double z, double dz) {
return (6.*z-6.*pow(z,2.0))/dz;
}
double f4z(double z) {
return 3.*pow(z,2.0)-2.*z;
}
double f1zz(double z,double dz) {
return (-6.+12.*z)/pow(dz,2.0);
}
double f2zz(double z,double dz) {
return (-4.+6.*z)/dz;
}
double f3zz(double z,double dz) {
return (6.-12.*z)/pow(dz,2.0);
}
double f4zz(double z,double dz) {
return (6.*z-2.)/dz;
}

void fPHI(double *solution, double gp1,double gp2,double dx, double dy) {
	solution[0]=f1(gp1)*f1(gp2);
	solution[1]=f2(gp1,dx)*f1(gp2);
	solution[2]=f1(gp1)*f2(gp2,dy);
	solution[3]=f2(gp1,dx)*f2(gp2,dy);
	solution[4]=f1(gp1)*f3(gp2);
	solution[5]=f2(gp1,dx)*f3(gp2);
	solution[6]=f1(gp1)*f4(gp2,dy);
	solution[7]=f2(gp1,dx)*f4(gp2,dy);
	solution[8]=f3(gp1)*f1(gp2);
	solution[9]=f4(gp1,dx)*f1(gp2);
	solution[10]=f3(gp1)*f2(gp2,dy);
	solution[11]=f4(gp1,dx)*f2(gp2,dy);
	solution[12]=f3(gp1)*f3(gp2);
	solution[13]=f4(gp1,dx)*f3(gp2);
	solution[14]=f3(gp1)*f4(gp2,dy);
	solution[15]=f4(gp1,dx)*f4(gp2,dy);
}

void fPHIX(double *solution, double gp1,double gp2,double dx,double dy) {
	solution[0]=f1z(gp1,dx)*f1(gp2);
	solution[1]=f2z(gp1)*f1(gp2);
	solution[2]=f1z(gp1,dx)*f2(gp2,dy);
	solution[3]=f2z(gp1)*f2(gp2,dy);
	solution[4]=f1z(gp1,dx)*f3(gp2);
	solution[5]=f2z(gp1)*f3(gp2);
	solution[6]=f1z(gp1,dx)*f4(gp2,dy);
	solution[7]=f2z(gp1)*f4(gp2,dy);
	solution[8]=f3z(gp1,dx)*f1(gp2);
	solution[9]=f4z(gp1)*f1(gp2);
	solution[10]=f3z(gp1,dx)*f2(gp2,dy);
	solution[11]=f4z(gp1)*f2(gp2,dy);
	solution[12]=f3z(gp1,dx)*f3(gp2);
	solution[13]=f4z(gp1)*f3(gp2);
	solution[14]=f3z(gp1,dx)*f4(gp2,dy);
	solution[15]=f4z(gp1)*f4(gp2,dy);
}

void fPHIY(double *solution, double gp1,double gp2,double dx,double dy) {
	solution[0]=f1(gp1)*f1z(gp2,dy);
	solution[1]=f2(gp1,dx)*f1z(gp2,dy);
	solution[2]=f1(gp1)*f2z(gp2);
	solution[3]=f2(gp1,dx)*f2z(gp2);
	solution[4] =f1(gp1)*f3z(gp2,dy);
	solution[5]=f2(gp1,dx)*f3z(gp2,dy);
	solution[6]=f1(gp1)*f4z(gp2);
	solution[7]=f2(gp1,dx)*f4z(gp2);
	solution[8]=f3(gp1)*f1z(gp2,dy);
	solution[9]=f4(gp1,dx)*f1z(gp2,dy);
	solution[10]=f3(gp1)*f2z(gp2);
	solution[11]=f4(gp1,dx)*f2z(gp2);
	solution[12]=f3(gp1)*f3z(gp2,dy);
	solution[13]=f4(gp1,dx)*f3z(gp2,dy);
	solution[14]=f3(gp1)*f4z(gp2);
	solution[15]=f4(gp1,dx)*f4z(gp2);
}

void fPHIXY(double *solution, double gp1,double gp2,double dx,double dy) {
	solution[0]=f1z(gp1,dx)*f1z(gp2,dy);
	solution[1]=f2z(gp1)*f1z(gp2,dy);
	solution[2]=f1z(gp1,dx)*f2z(gp2);
	solution[3]=f2z(gp1)*f2z(gp2);
	solution[4]=f1z(gp1,dx)*f3z(gp2,dy);
	solution[5]=f2z(gp1)*f3z(gp2,dy);
	solution[6]=f1z(gp1,dx)*f4z(gp2);
	solution[7]=f2z(gp1)*f4z(gp2);
	solution[8]=f3z(gp1,dx)*f1z(gp2,dy);
	solution[9]=f4z(gp1)*f1z(gp2,dy);
	solution[10]= f3z(gp1,dx)*f2z(gp2);
	solution[11]=f4z(gp1)*f2z(gp2);
	solution[12]=f3z(gp1,dx)*f3z(gp2,dy);
	solution[13]=f4z(gp1)*f3z(gp2,dy);
	solution[14]=f3z(gp1,dx)*f4z(gp2);
	solution[15]=f4z(gp1)*f4z(gp2);
}

void fPHIXX(double *solution, double gp1,double gp2,double dx,double dy) {
	solution[0]=f1zz(gp1,dx)*f1(gp2);
	solution[1]=f2zz(gp1,dx)*f1(gp2);
	solution[2]=f1zz(gp1,dx)*f2(gp2,dy);
	solution[3]=f2zz(gp1,dx)*f2(gp2,dy);
	solution[4]=f1zz(gp1,dx)*f3(gp2);
	solution[5]=f2zz(gp1,dx)*f3(gp2);
	solution[6]=f1zz(gp1,dx)*f4(gp2,dy);
	solution[7]=f2zz(gp1,dx)*f4(gp2,dy);
	solution[8]=f3zz(gp1,dx)*f1(gp2);
	solution[9]=f4zz(gp1,dx)*f1(gp2);
	solution[10]=f3zz(gp1,dx)*f2(gp2,dy);
	solution[11]=f4zz(gp1,dx)*f2(gp2,dy);
	solution[12]=f3zz(gp1,dx)*f3(gp2);
	solution[13]=f4zz(gp1,dx)*f3(gp2);
	solution[14]=f3zz(gp1,dx)*f4(gp2,dy);
	solution[15]=f4zz(gp1,dx)*f4(gp2,dy);
}

void fPHIYY(double *solution, double gp1,double gp2,double dx,double dy) {
	solution[0]=f1(gp1)*f1zz(gp2,dy);
	solution[1]=f2(gp1,dx)*f1zz(gp2,dy);
	solution[2]=f1(gp1)*f2zz(gp2,dy);
	solution[3]=f2(gp1,dx)*f2zz(gp2,dy);
	solution[4]=f1(gp1)*f3zz(gp2,dy);
	solution[5]=f2(gp1,dx)*f3zz(gp2,dy);
	solution[6]=f1(gp1)*f4zz(gp2,dy);
	solution[7]=f2(gp1,dx)*f4zz(gp2,dy);
	solution[8]=f3(gp1)*f1zz(gp2,dy);
	solution[9]=f4(gp1,dx)*f1zz(gp2,dy);
	solution[10]=f3(gp1)*f2zz(gp2,dy);
	solution[11]=f4(gp1,dx)*f2zz(gp2,dy);
	solution[12]=f3(gp1)*f3zz(gp2,dy);
	solution[13]=f4(gp1,dx)*f3zz(gp2,dy);
	solution[14]=f3(gp1)*f4zz(gp2,dy);
	solution[15]=f4(gp1,dx)*f4zz(gp2,dy);
}
