function coef_n2=n2_characterization(distribution_type,ne,n2,...
    spec,nex,ney,xlen,ylen)

% distribution_type==1 Polynomical in x-dir
% distribution_type==2 Polynomical in y-dir

gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
dx=xlen/nex;
dy=ylen/ney;

coef_n2=zeros(ne,3,3,3);

if distribution_type==1
    
    if length(spec)~=2 || spec(2)>=xlen/2
        error('spec=[degree, width of flat area]')
    end
    gp=gp*dx;
    x1=spec(2);
    x2=xlen-x1;
    degree_plus=spec(1)+1;
    mx=mx_polynomial_order_constraint(spec(1),x1,x2);
    b=zeros(degree_plus,1);
    b(1)=n2(1);
    b(2)=n2(2);
    coef=mx\b;
    multiplier1=zeros(1,degree_plus);
    multiplier2=zeros(1,degree_plus);
    for e=1:1:ne
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            if xPos<=x1
                n2_abs=n2(1);
                dn2_dx=0;
            elseif xPos<=x2 
                for i=1:1:degree_plus
                    multiplier1(i)=xPos^(degree_plus-i);
                end
                for i=1:1:spec(1)
                    multiplier2(i)=(degree_plus-i)*xPos^(spec(1)-i);
                end
                n2_abs=dot(coef,multiplier1);
                dn2_dx=dot(coef,multiplier2);
            elseif xPos<=xlen
                n2_abs=n2(2);
                dn2_dx=0;
            end
            for iy=1:1:3
                coef_n2(e,ix,iy,1)=n2_abs;
                coef_n2(e,ix,iy,2)=dn2_dx;
            end
        end
    end
    
elseif distribution_type==2
    
    if length(spec)~=2 || spec(2)>=ylen/2
        error('spec=[degree, width of flat area]')
    end
    gp=gp*dy;
    y1=spec(2);
    y2=ylen-y1;
    degree_plus=spec(1)+1;
    mx=mx_polynomial_order_constraint(spec(1),y1,y2);
    b=zeros(degree_plus,1);
    b(1)=n2(1);
    b(2)=n2(2);
    coef=mx\b;
    multiplier1=zeros(1,degree_plus);
    multiplier2=zeros(1,degree_plus);
    for e=1:1:ne
        [~,iny]=inxiny_elemental(e,ney);
        yPos0=(iny-1)*dy;
        for iy=1:1:3
            yPos=yPos0+gp(iy);
            if yPos<=y1
                n2_abs=n2(1);
                dn2_dy=0;
            elseif yPos<=y2 
                for i=1:1:degree_plus
                    multiplier1(i)=yPos^(degree_plus-i);
                end
                for i=1:1:spec(1)
                    multiplier2(i)=(degree_plus-i)*yPos^(spec(1)-i);
                end
                n2_abs=dot(coef,multiplier1);
                dn2_dy=dot(coef,multiplier2);
            elseif yPos<=ylen
                n2_abs=n2(2);
                dn2_dy=0;
            end
            for ix=1:1:3
                coef_n2(e,ix,iy,1)=n2_abs;
                coef_n2(e,ix,iy,3)=dn2_dy;
            end
        end
    end
    
else
    error('distribution_type specified is not available')

end

end

% function mx=mx_polynomial_order_constraint(degree,x1,x2)
% if mod(degree,2)~=1
%     error('Degree must be an odd number')
% end
% degree_plus=degree+1;
% mx=zeros(degree_plus,degree_plus);
% for order=1:1:degree_plus/2
%     for column=1:1:degree_plus
%         deg=degree_plus-column-order+1;
%         if order==1
%             coef=1;
%         else
%             coef=factorial_sp(degree_plus-column,order-1);
%         end
%         pre_index=2*order-2;
%         if deg<0
%             mx(pre_index+1,column)=0;
%         else
%             mx(pre_index+1,column)=coef*x1^deg;
%         end
%         
%         mx(pre_index+2,column)=coef*x2^deg;
%     end
% end
% 
% end
% 
% function sol=factorial_sp(val,n)
% sol=val;
% for i=1:1:n-1
%     sol=sol*(val-i);
% end
% end