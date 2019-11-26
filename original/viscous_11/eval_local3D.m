function [con,cono,conx,cony,conz,conxx,conyy,conzz]=eval_local3D(c,co,phi,phix,phiy,phiz,phixx,phiyy,phizz,gbasis)

con=0.0;
cono=0.0;
conx=0.0;
cony=0.0;
conz=0.0;
conxx=0.0;
conyy=0.0;
conzz=0.0;

for i=1:64
    con=con+c(gbasis(i))*phi(i);
    cono=cono+co(gbasis(i))*phi(i);
    conx=conx+c(gbasis(i))*phix(i);
    cony=cony+c(gbasis(i))*phiy(i);
    conz=conz+c(gbasis(i))*phiz(i);
    conxx=conxx+c(gbasis(i))*phixx(i);
    conyy=conyy+c(gbasis(i))*phiyy(i);
    conzz=conzz+c(gbasis(i))*phizz(i);
end