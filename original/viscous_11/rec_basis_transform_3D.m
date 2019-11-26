function [phi,phix,phiy,phiz,phixx,phiyy,phizz]=rec_basis_transform_3D(dx,dy,dz,phi,phix,phiy,phiz,phixx,phiyy,phizz)

type=...
    [0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    1 0 1
    0 1 1
    1 1 1];

k=0;
for iorientation=1:1:8
    for itype=1:1:8
        k=k+1;
        if type(itype,1)==0
            phix(k)=phix(k)/dx;
            phixx(k)=phixx(k)/(dx^2);
        else
            phi(k)=phi(k)*dx;
            phixx(k)=phixx(k)/dx;
        end
        if type(itype,2)==0
            phiy(k)=phiy(k)/dy;
            phiyy(k)=phiyy(k)/(dy^2);
        else
            phi(k)=phi(k)*dy;
            phiyy(k)=phiyy(k)/dy;
        end
        if type(itype,3)==0
            phiz(k)=phiz(k)/dz;
            phizz(k)=phizz(k)/(dz^2);
        else
            phi(k)=phi(k)*dz;
            phizz(k)=phizz(k)/dz;
        end
    end
end

end