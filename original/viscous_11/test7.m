function test7

clear
clc

nex=5;
ne=nex;
w=[0.27778 0.4444 0.27778];
gp=[0.1127 0.5 0.8873];
nx=nex+1;
n=nx;
nfour=n*4;
c=ones(1,nfour)*2;
c=bc_septic(c);
% c=co;
x=linspace(0,1,nx);
for rahmat=1:1:1   
    sf=zeros(1,nfour);
    sj=zeros(nfour,nfour);
    for e=1:1:ne
        dx=x(e+1)-x(e);
        for w_=1:1:3
            phi=get_phi_septic(gp(w_),dx);
            phixxxx=get_phixxxx_septic(gp(w_),dx);
            yxxxx=0.0;
            for i=1:1:8
                yxxxx=yxxxx+c(e*4-4+i)*phixxxx(i);
            end
            for i=1:1:8 
                sf(e*4-4+i)=sf(e*4-4+i)-w(w_)*dx*(yxxxx-360*gp(w_)^2-48)*phi(i);
                for j=1:1:8
                    if w(w_)*dx*phi(i)*phi(j)==0
                        error('jhgjk')
                    end
                    sj(e*4-4+i,e*4-4+j)=sj(e*4-4+i,e*4-4+j)+w(w_)*dx*phi(i)*phixxxx(j);
                end
            end
        end
    end
    sf=bcsf_septic(sf);
    sj=bcsj_septic(sj,nfour);
    det(sj)
    if abs(det(sj))<=1
        sj
        ddd=diag(sj)
        abs(diag(sj))
        [val,ii]=min(abs(diag(sj)))
        val
        ii
        ddd(ii)
        
        det(sj)
        error('Singular')
    end
    c_=sj\sf';
    c=c+c_';
    err = sqrt(sum(c_.^2.0))
end
sj
c__=zeros(1,n);
for i=1:1:n
    c__(i)=c(4*i-3);
end
% plot(x,c__)
n=100;
x=linspace(0,1,n);
plot(x,x.^6+2*x.^4-5*x.^2+6*ones(1,n))
grid on

end

function phi=get_phi_septic(val,dx)
phi=zeros(1,8);
k=0;
for orientation=[0 1]
    for type=0:1:3
        k=k+1;
        phi(k)=septic_Hermite_basis(val,orientation,type,0)*(dx^(type));
    end
end
end

function phixxxx=get_phixxxx_septic(val,dx)
phixxxx=zeros(1,8);
k=0;
for orientation=[0 1]
    for type=0:1:3
        k=k+1;
        phixxxx(k)=septic_Hermite_basis(val,orientation,type,4)*(dx^(type))/(dx^4);
    end
end
end

function c=bc_septic(c)
c(1)=6;
c(2)=0;
c(3)=-10;
c(4)=0;
end

function sf=bcsf_septic(sf)
sf(1)=0;
sf(2)=0;
sf(3)=0;
sf(4)=0;
end

function sj=bcsj_septic(sj,nfour)
for i=1:1:nfour
    sj(1,i)=0;
    sj(2,i)=0;
    sj(3,i)=0;
    sj(4,i)=0;
end
sj(1,1)=1;
sj(2,2)=1;
sj(3,3)=1;
sj(4,4)=1;
end
