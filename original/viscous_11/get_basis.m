function [phi,phix,phiy,phiz,phixx,phiyy,phizz]=get_basis(alpha,beta,gamma)

%alpha, beta, and gamma are in local coordinates

phi=zeros(1,64);
phix=zeros(1,64);
phiy=zeros(1,64);
phiz=zeros(1,64);
phixx=zeros(1,64);
phiyy=zeros(1,64);
phizz=zeros(1,64);

orientation=...
    [0 0 0
    0 1 0
    1 0 0
    1 1 0
    0 0 1
    0 1 1
    1 0 1
    1 1 1];

type=...
    [0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    1 0 1
    0 1 1
    1 1 1];

i=0;
for iorientation=1:1:8
    for itype=1:1:8
        i=i+1;
        
        phi(i)=basis(alpha,orientation(iorientation,1),type(itype,1),0)*...
            basis(beta,orientation(iorientation,2),type(itype,2),0)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),0);
        
        phix(i)=basis(alpha,orientation(iorientation,1),type(itype,1),1)*...
            basis(beta,orientation(iorientation,2),type(itype,2),0)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),0);
        
        phiy(i)=basis(alpha,orientation(iorientation,1),type(itype,1),0)*...
            basis(beta,orientation(iorientation,2),type(itype,2),1)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),0);
        
        phiz(i)=basis(alpha,orientation(iorientation,1),type(itype,1),0)*...
            basis(beta,orientation(iorientation,2),type(itype,2),0)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),1);
        
        phixx(i)=basis(alpha,orientation(iorientation,1),type(itype,1),2)*...
            basis(beta,orientation(iorientation,2),type(itype,2),0)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),0);
        
        phiyy(i)=basis(alpha,orientation(iorientation,1),type(itype,1),0)*...
            basis(beta,orientation(iorientation,2),type(itype,2),2)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),0);
        
        phizz(i)=basis(alpha,orientation(iorientation,1),type(itype,1),0)*...
            basis(beta,orientation(iorientation,2),type(itype,2),0)*...
            basis(gamma,orientation(iorientation,3),type(itype,3),2);
    end
end

end