function illustrate_mass_flux(c,nx,ny,n1,n2,T,diff,conc_,coef_T)
%ONLY FOR LINEAR TEMPERATURE!
%---------------------
slope=T(2)-T(1);
dx=1/(nx-1);
temp=T(1)-slope*dx;
%----------------------

% size(coef_T)
% error('fdafdasg')


c_nodal=extract_nodal_weights_2D(c,nx,ny);
[DX,DY]=extract_gradient_2D(c,nx,ny);


[DXX,~]=gradient(DX);
[~,DYY]=gradient(DY);


%----------------------------------------
w=[5/18 4/9 5/18];
% T=T(1);
% chi=0.5-(1-1/T);
nex=nx-1;
ney=ny-1;
ne=nex*ney;
energy=zeros(ney,nex);
e=0;
dxdy=1/(nex*ney);
dx=1/nex;
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2]/nex;
for xth=1:1:nex
    for yth=1:1:ney
        e=e+1;
%         [inx,~]=inxiny_elemental(e,ney);
%         pos0=(inx-1)*dx;
        for ix=1:1:3
%             T_=T(1)+slope*(pos0+gp(ix));
%             chi=0.5-(1-1/T_);
            for iy=1:1:3
                T_=coef_T(e,ix,iy,1);
                chi=0.5-(1-1/T_);
                c2=1-conc_(ix,iy,1,e);
                first=diff*T_*((log(conc_(ix,iy,1,e))+1)/n1-(log(c2)+1)/n2+chi*(1-2*conc_(ix,iy,1,e))/n1);
                
%                 first=T_;
                
                second=conc_(ix,iy,4,e)+conc_(ix,iy,5,e);
                
%                 second=0;
                
                energy(yth,xth)=energy(yth,xth)+w(ix)*w(iy)...
                    *(first-second);
                
                first/second;
            end
        end
        energy(yth,xth)=energy(yth,xth)*dxdy;
    end
end

subplot(2,2,1)
surf(c_nodal)
subplot(2,2,2)
surf(energy)

[DX,DY]=gradient(energy);

subplot(2,2,3)
surf(DX)
subplot(2,2,4)
surf(DY)

return

T=T(1);
chi=0.5-(1-1/T);
P=zeros(ny,nx);
P1=P;P2=P;
for xth=1:1:nx
    for yth=1:1:ny
        c2=1-c_nodal(yth,xth);
        P1(yth,xth)=diff*T...
            *(c_nodal(yth,xth)*log(c_nodal(yth,xth))/n1...
            +c2*log(c2)/n2+c_nodal(yth,xth)*c2*chi/n1);
        P2(yth,xth)=...
            0.5...
            *(DX(yth,xth)^2+DY(yth,xth)^2);
        P(yth,xth)=P1(yth,xth)+P2(yth,xth);
    end
end



[DX_,DY_]=gradient(c_nodal);

sum(sum(abs(DX)))
sum(sum(abs(DY)))

sum(sum(abs(DX_)))
sum(sum(abs(DY_)))

DX_diff=DX-DY_;
DY_diff=DY-DY_;

% size(DX_diff)
% error('dfaas')

sum(sum(abs(DX_diff)))
sum(sum(abs(DY_diff)))

subplot(2,3,1)
surf(c_nodal)
subplot(2,3,2)
surf(DX)
subplot(2,3,3)
surf(DY)
subplot(2,3,4)
surf(P1)
subplot(2,3,5)
surf(P2)
subplot(2,3,6)
surf(P)


return

error('fafaf')
%-------------------------------------------

laplacian=DXX+DYY;

chunk=zeros(ny,nx);
for xth=1:1:nx
    %--------------
    temp=temp+slope*dx;
    chi=0.5-(1-1/temp);
    %------------
    for yth=1:1:ny
        con=c_nodal(yth,xth);
        con2=1-con;
        dfdc=(log(con)+1)/n1-(log(con2)+1)/n2+chi*(1-2*con)/n1;
        chunk(yth,xth)=diff*temp*dfdc-laplacian(yth,xth);
    end
end

[xFlux,yFlux]=gradient(chunk);

subplot(2,2,1)
surf(c_nodal);
subplot(2,2,2)
surf(laplacian)
subplot(2,2,3)
surf(DX);
subplot(2,2,4)
surf(DY)


end