function [coef_T,two_chi_n1]=generate_T_dist_3d(...
    type,spec,T1,T2,ne,nex,ney,nez,xlen,ylen,zlen,entropy,n1)


if type==1
    [coef_T,two_chi_n1]=generate_linear_T_grad_3d(...
        T1,T2,ne,nex,ney,nez,zlen,entropy,n1);
elseif type==2
    [coef_T,two_chi_n1]=generate_polynomial_T_grad_3d(...
        spec,T1,T2,ne,nex,ney,nez,zlen,entropy,n1);
elseif type==3
    [coef_T,two_chi_n1]=generate_sinusoidal_T_3d(...
        ne,T1,T2,xlen,ylen,nex,ney,nez,entropy,n1);
elseif type==4
    [coef_T,two_chi_n1]=generate_donut_T_3d(...
        spec,ne,T1,T2,xlen,ylen,zlen,nex,ney,nez,entropy,n1);
else
    error('Specified type does not exist')
end

end

function [coef_T,two_chi_n1]=generate_donut_T_3d(...
    spec,ne,T1,T2,xlen,ylen,zlen,nex,ney,nez,entropy,n1)

coef_T=zeros(ne,3,3,3,4);
two_chi_n1=zeros(ne,3,3,3);
% A=(T2-T1)/4;
wx=2*pi/xlen;
wy=2*pi/ylen;
nexney=nex*ney;
dx=xlen/nex;
dy=ylen/ney;
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
gpx=gp*dx;
gpy=gp*dy;
% z1=nexney*(nez-1);
iny=0;


%----------------------------------

wz=2*pi*spec/zlen;
A=(T2-T1)/8;
inz=0;
dz=zlen/nez;
gpz=gp*dz;
%----------------------------------

for e=1:1:ne
    
    if mod(e,nexney)==1
        inz=inz+1;
        zPos0=(inz-1)*dz;
    end

    if mod(e,nex)==1
        inx=0;
        iny=iny+1;
    end
    inx=inx+1;
    xPos0=(inx-1)*dx;
    yPos0=(iny-1)*dy;
    for ix=1:1:3
        xPos=xPos0+gpx(ix);
        for iy=1:1:3
            yPos=yPos0+gpy(iy);
            
            for iz=1:1:3
                
                zPos=zPos0+gpz(iz);
            
                temp=A*((1-cos(wx*xPos))*(1-cos(wy*yPos))*(1-cos(wz*zPos)))+T1;


                dT_dx=A*((wx*sin(wx*xPos))*(1-cos(wy*yPos))*(1-cos(wz*zPos)));
                dT_dy=A*((1-cos(wx*xPos))*(wy*sin(wy*yPos))*(1-cos(wz*zPos)));
                dT_dz=A*((1-cos(wx*xPos))*(1-cos(wy*yPos))*(wz*sin(wz*zPos)));

                two_chi_n1_=2*dimensionless_chi(temp,entropy)/n1;
                
                coef_T(e,ix,iy,iz,1)=temp;
                coef_T(e,ix,iy,iz,2)=dT_dx;
                coef_T(e,ix,iy,iz,3)=dT_dy;
                coef_T(e,ix,iy,iz,4)=dT_dz;
                two_chi_n1(e,ix,iy,iz)=two_chi_n1_;

    %             for iz=1:1:3
    %                 coef_T(e,ix,iy,iz,1)=temp;
    %                 coef_T(e,ix,iy,iz,2)=dT_dx;
    %                 coef_T(e,ix,iy,iz,3)=dT_dy;
    %                 two_chi_n1(e,ix,iy,iz)=two_chi_n1_;
    %             end
            end
        end
    end
    
%     for e_=e:nexney:e+z1
%         for ix=1:1:3
%             for iy=1:1:3
%                 for iz=1:1:3
%                     coef_T(e_,ix,iy,iz,1)=coef_T(e,ix,iy,iz,1);
%                     coef_T(e_,ix,iy,iz,2)=coef_T(e,ix,iy,iz,2);
%                     coef_T(e_,ix,iy,iz,3)=coef_T(e,ix,iy,iz,3);
%                     two_chi_n1(e_,ix,iy,iz)=two_chi_n1(e,ix,iy,iz);
%                 end
%             end
%         end
%     end
    
end

end

function [coef_T,two_chi_n1]=generate_sinusoidal_T_3d(...
    ne,T1,T2,xlen,ylen,nex,ney,nez,entropy,n1)

coef_T=zeros(ne,3,3,3,4);
two_chi_n1=zeros(ne,3,3,3);
A=(T2-T1)/4;
wx=2*pi/xlen;
wy=2*pi/ylen;
nexney=nex*ney;
dx=xlen/nex;
dy=ylen/ney;
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
gpx=gp*dx;
gpy=gp*dy;
z1=nexney*(nez-1);
iny=0;
for e=1:1:nexney

    if mod(e,nex)==1
        inx=0;
        iny=iny+1;
    end
    inx=inx+1;
    xPos0=(inx-1)*dx;
    yPos0=(iny-1)*dy;
    for ix=1:1:3
        xPos=xPos0+gpx(ix);
        for iy=1:1:3
            yPos=yPos0+gpy(iy);
            
            temp=A*((1-cos(wx*xPos))*(1-cos(wy*yPos)))+T1;
            dT_dx=A*((wx*sin(wx*xPos))*(1-cos(wy*yPos)));
            dT_dy=A*((1-cos(wx*xPos))*(wy*sin(wy*yPos)));
            
            two_chi_n1_=2*dimensionless_chi(temp,entropy)/n1;
            
            for iz=1:1:3
                coef_T(e,ix,iy,iz,1)=temp;
                coef_T(e,ix,iy,iz,2)=dT_dx;
                coef_T(e,ix,iy,iz,3)=dT_dy;
                two_chi_n1(e,ix,iy,iz)=two_chi_n1_;
            end
            
%             temperature(e,ix,iy,1)=A*((1-cos(c*xPos))*(1-cos(c*yPos)))+T(1);
%             temperature(e,ix,iy,2)=A*((c*sin(c*xPos))*(1-cos(c*yPos)));
%             temperature(e,ix,iy,3)=A*((1-cos(c*xPos))*(c*sin(c*yPos)));
        end
    end
    
    for e_=e:nexney:e+z1
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    coef_T(e_,ix,iy,iz,1)=coef_T(e,ix,iy,iz,1);
                    coef_T(e_,ix,iy,iz,2)=coef_T(e,ix,iy,iz,2);
                    coef_T(e_,ix,iy,iz,3)=coef_T(e,ix,iy,iz,3);
                    two_chi_n1(e_,ix,iy,iz)=two_chi_n1(e,ix,iy,iz);
                end
            end
        end
    end
    
end

end

function [coef_T,two_chi_n1]=generate_polynomial_T_grad_3d(...
    spec,T1,T2,ne,nex,ney,nez,zlen,entropy,n1)

if length(spec)~=2 || spec(2)>=zlen/2
    error('T_distribution_spec=[degree, width of flat area]')
end

z1=spec(2);
z2=zlen-z1;
degree_plus=spec(1)+1;
mx=mx_polynomial_order_constraint(spec(1),z1,z2);
b=zeros(degree_plus,1);
b(1)=T1;
b(2)=T2;
coef=mx\b;

coef_T=zeros(ne,3,3,3,4);
two_chi_n1=zeros(ne,3,3,3);

zth_=-1;
nexney=nex*ney;
dz=zlen/nez;

dz_minor=dz*[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
slope=(T2-T1)/zlen;

multiplier1=zeros(1,degree_plus);
multiplier2=zeros(1,degree_plus);

for e=1:1:ne
    if mod(e,nexney)==1
        zth_=zth_+1;
        pre_z=zth_*dz;
    end
    
    for zth=1:1:3
        
        zPos=pre_z+dz_minor(zth);
        
        if zPos<=z1
            T_abs=T1;
            dT_dz=0;
        
        elseif zPos<z2
%             T=T1+slope*(pre_z+dz_minor(zth));
%             two_chi_n1_=2*dimensionless_chi(T,entropy)/n1;
%             for yth=1:1:3
%                 for xth=1:1:3
%                     coef_T(e,xth,yth,zth,1)=T;
%                     coef_T(e,xth,yth,zth,4)=slope;
%                     two_chi_n1(e,xth,yth,zth)=two_chi_n1_;
%                 end
%             end
            for i=1:1:degree_plus
                multiplier1(i)=zPos^(degree_plus-i);
            end
            T_abs=dot(coef,multiplier1);
            for i=1:1:spec(1)     
                multiplier2(i)=(degree_plus-i)*zPos^(spec(1)-i);
            end
            dT_dz=dot(coef,multiplier2);
            
        elseif zPos<=zlen
            T_abs=T2;
            dT_dz=0;
        end
        
        two_chi_n1_=2*dimensionless_chi(T_abs,entropy)/n1;
        for yth=1:1:3
            for xth=1:1:3
                coef_T(e,xth,yth,zth,1)=T_abs;
                coef_T(e,xth,yth,zth,4)=dT_dz;
                two_chi_n1(e,xth,yth,zth)=two_chi_n1_;
            end
        end
        
    end
    
end

end

function [coef_T,two_chi_n1]=generate_linear_T_grad_3d(...
    T1,T2,ne,nex,ney,nez,zlen,entropy,n1)

coef_T=zeros(ne,3,3,3,4);
two_chi_n1=zeros(ne,3,3,3);

zth_=-1;
nexney=nex*ney;
dz=zlen/nez;

dz_minor=dz*[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
slope=(T2-T1)/zlen;

for e=1:1:ne
    if mod(e,nexney)==1
        zth_=zth_+1;
        pre_z=zth_*dz;
    end
    
    for zth=1:1:3
        T=T1+slope*(pre_z+dz_minor(zth));
        two_chi_n1_=2*dimensionless_chi(T,entropy)/n1;
        for yth=1:1:3
            for xth=1:1:3
                coef_T(e,xth,yth,zth,1)=T;
                coef_T(e,xth,yth,zth,4)=slope;
                two_chi_n1(e,xth,yth,zth)=two_chi_n1_;
            end
        end
    end
    
end

end