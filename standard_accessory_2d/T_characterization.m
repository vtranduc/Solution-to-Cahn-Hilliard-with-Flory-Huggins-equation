function temperature=T_characterization(distribution_type,ne,T,...
    T_distribution_spec,nex,ney,xlen,ylen,...
    diff,n1,n2,ci_ave,chi_a,chi_b,chi_type)
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
dx=xlen/nex;
gp=gp*dx;
temperature=zeros(ne,3,3,3);
if distribution_type==1
    % Polynomial 3 degrees in x-dir
    if length(T_distribution_spec)~=1 || T_distribution_spec<0 || T_distribution_spec>=0.5
        error('T_distribution_spec must be larger or equal to 0 and smaller than 0.5!')
    end
    pos2=xlen-T_distribution_spec;
    divisor=[T(1);T(2);0;0];
    M=[T_distribution_spec^3 T_distribution_spec^2 T_distribution_spec 1
        pos2^3 pos2^2 pos2 1
        3*T_distribution_spec^2 2*T_distribution_spec 1 0
        3*pos2^2 2*pos2 1 0];
    coef=M\divisor;
    
%     coef
    
    for e=1:1:ne
        [inx,~]=inxiny_elemental(e,ney);
        ifirst=(inx-1)*dx;
        for einx=1:1:3
            pos=ifirst+gp(einx);
            if pos<T_distribution_spec
                for einy=1:1:3
                    temperature(e,einx,einy,1)=T(1);
                end
            elseif pos>pos2
                for einy=1:1:3
                    temperature(e,einx,einy,1)=T(2);
                end
            else
                temperature(e,einx,1,1)=dot(coef,[pos^3 pos^2 pos 1]);
                temperature(e,einx,1,2)=3*coef(1)*pos^2+2*coef(2)*pos+coef(3);
%                 temperature(e,einx,1,4)=6*coef(1)*pos+2*coef(2);
                for einy=2:3
                    temperature(e,einx,einy,1)=temperature(e,einx,1,1);
                    temperature(e,einx,einy,2)=temperature(e,einx,1,2);
%                     temperature(e,einx,einy,4)=temperature(e,einx,1,4);
                end
            end
        end
    end
    %==================TEST TYPES============================================
elseif distribution_type==0
    % Linear in x dir
    slope=(T(2)-T(1))/xlen;
    for e=1:1:ne
        [xth,~]=inxiny_elemental(e,ney);
        ifirst=(xth-1)*dx;
        for ix=1:1:3
            pos=ifirst+gp(ix);
            temperature(e,ix,1,1)=slope*pos+T(1);
            temperature(e,ix,1,2)=slope;
            for iy=2:1:3
                temperature(e,ix,iy,1)=temperature(e,ix,1,1);
                temperature(e,ix,iy,2)=slope;
            end
        end
    end
    %==================================================================
elseif distribution_type==2
    %Polynomial of 5 degrees in x dir
    pos2=xlen-T_distribution_spec;
    M=...
        [T_distribution_spec^5 T_distribution_spec^4 T_distribution_spec^3 T_distribution_spec^2 T_distribution_spec 1
        pos2^5 pos2^4 pos2^3 pos2^2 pos2 1
        5*T_distribution_spec^4 4*T_distribution_spec^3 3*T_distribution_spec^2 2*T_distribution_spec 1 0
        5*pos2^4 4*pos2^3 3*pos2^2 2*pos2 1 0
        20*T_distribution_spec^3 12*T_distribution_spec^2 6*T_distribution_spec 2 0 0
        20*pos2^3 12*pos2^2 6*pos2 2 0 0];
    coef=M\[T(1);T(2);0;0;0;0];  
    for e=1:1:ne
        [inx,~]=inxiny_elemental(e,ney);
        pos0=(inx-1)*dx;
        for ix=1:1:3
            pos=pos0+gp(ix);
            if pos<T_distribution_spec
                for iy=1:1:3
                    temperature(e,ix,iy,1)=T(1);
                end
            elseif pos>pos2
                for iy=1:1:3
                    temperature(e,ix,iy,1)=T(2);
                end
            else
                temperature(e,ix,1,1)=...
                    coef(1)*pos^5+coef(2)*pos^4+coef(3)*pos^3+coef(4)*pos^2+coef(5)*pos+coef(6);
                temperature(e,ix,1,2)=...
                    5*coef(1)*pos^4+4*coef(2)*pos^3+3*coef(3)*pos^2+2*coef(4)*pos+coef(5);
%                 temperature(e,ix,1,4)=...
%                     20*coef(1)*pos^3+12*coef(2)*pos^2+6*coef(3)*pos+2*coef(4);
                for iy=2:3
                    temperature(e,ix,iy,1)=temperature(e,ix,1,1);
                    temperature(e,ix,iy,2)=temperature(e,ix,1,2);
%                     temperature(e,ix,iy,4)=temperature(e,ix,1,4);
                end
            end
        end
    end
elseif distribution_type==3
    %Cosine wave centered in the middle
    
    gpy=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2]*dy;
    
    A=(T(2)-T(1))/4;
    
    wx=2*pi/xlen;
    wy=2*pi/ylen;
    
    dy=ylen/ney;
    for e=1:1:ne
        [inx,iny]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        yPos0=(iny-1)*dy;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            for iy=1:1:3
                yPos=yPos0+gpy(iy);
                temperature(e,ix,iy,1)=A*((1-cos(wx*xPos))*(1-cos(wy*yPos)))+T(1);
                temperature(e,ix,iy,2)=A*((wx*sin(wx*xPos))*(1-cos(wy*yPos)));
                temperature(e,ix,iy,3)=A*((1-cos(wx*xPos))*(wy*sin(wy*yPos)));
            end
        end
    end
elseif distribution_type==4
    %Radial concentric cicle
    if length(T_distribution_spec)~=2 || T_distribution_spec(1)>=T_distribution_spec(2) || T_distribution_spec(2)>0.5 || T_distribution_spec(1)<0
        error('T_distribution_spec must have 2 numbers, where T_distribution_spec(1)<T_distribution_spec(2) && T_distribution_spec(2)<=0.5 && T_distribution_spec(1)>=0')
    end
    A=(T(2)-T(1))/2;
    dy=ylen/ney;
    z=pi/(T_distribution_spec(2)-T_distribution_spec(1));
    Az=-A*z;
    for e=1:1:ne
        [inx,iny]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        yPos0=(iny-1)*dy;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            for iy=1:1:3
                yPos=yPos0+gp(iy);
                r=sqrt((xPos-0.5)^2+(yPos-0.5)^2);
                if r>T_distribution_spec(1)
                    if r<T_distribution_spec(2)
                        dist=T_distribution_spec(2)-r;
                        temperature(e,ix,iy,1)=A*(1-cos(z*dist))+T(1);
                        temperature(e,ix,iy,2)=Az*sin(z*dist)*(xPos-0.5)/r;
                        temperature(e,ix,iy,3)=Az*sin(z*dist)*(yPos-0.5)/r;
                    elseif r>=T_distribution_spec(2)
                        temperature(e,ix,iy,1)=T(1);
                    end
                elseif r<=T_distribution_spec(1)
                    temperature(e,ix,iy,1)=T(2);
                end
            end
        end
    end
    %----------------------------------------
elseif distribution_type==5
    
    %Looks so wrong
    
    if length(T_distribution_spec)~=2 || T_distribution_spec(2)>ylen || T_distribution_spec(2)<0
        error('T_distribution_spec must have 2 numbers, and T_distribution_spec(2) must be within boundary')
    end
    dy=ylen/ney;
    T_ave=mean(T);
    Ao=(T(2)-T_ave)/2;
    
%     slope=Ao/T_distribution_spec(2);
    
    
    
    z1=T_distribution_spec(1)*2*pi;
    z2=pi/T_distribution_spec(2);
    z3=-z2*Ao;
    
    
    
    for e=1:1:ne
        
        [inx,iny]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        yPos0=(iny-1)*dy;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            for iy=1:1:3
                yPos=yPos0+gp(iy);
                
                if yPos>=T_distribution_spec(2)
                    temperature(e,ix,iy,1)=T_ave;
                elseif yPos<T_distribution_spec(2)
                    
                    A=Ao*(cos(z2*yPos)+1);
                    sine_term=sin(z1*xPos);
                    temperature(e,ix,iy,1)=A*sine_term+T_ave;

                    temperature(e,ix,iy,2)=z1*A*cos(z1*xPos);
                    temperature(e,ix,iy,3)=z3*sin(z2*yPos)*sine_term;
                    
                end
            end
        end
        
    end
    %--------------------------------------------
elseif distribution_type==6
    
    % Temperature varies with diagonal line from (0,0) to (xlen,ylen)
    % Use projection to map
    % Derivatives are obtained by directional derivatives
    
    len_sq_diag=xlen^2+ylen^2;
    x_=xlen/len_sq_diag;
    y_=ylen/len_sq_diag;
    dy=ylen/ney;
    len_diag=sqrt(len_sq_diag);
    mT=(T(2)-T(1))/len_diag;
    
    slopes=[xlen/len_diag,ylen/len_diag;-xlen/len_diag,ylen/len_diag]\[mT;0];
    
    for e=1:1:ne
        
        [inx,iny]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        yPos0=(iny-1)*dy;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            for iy=1:1:3
                yPos=yPos0+gp(iy); 
                multiplier=xPos*x_+yPos*y_;
                
                r1=multiplier*xlen;
                r2=multiplier*ylen;
                r=sqrt(r1^2+r2^2);
                temperature(e,ix,iy,1)=mT*r+T(1);
                
                temperature(e,ix,iy,2)=slopes(1);
                temperature(e,ix,iy,3)=slopes(2);
                
            end
        end
    end
elseif distribution_type==7
    
    % Flat at the edge, linear in the middle, polynomial in between
    
    % T_distribution_spec=[Width of flat temperature, width of polynomial temperature]
    
    if length(T_distribution_spec)~=2 || T_distribution_spec(1)+T_distribution_spec(2)>xlen/2
        error('T_distribution_spec must be 3. Firstif the length of flat area, second is the area of the curve. The must not add up to more than half of xlen')
    end
    
    x1=T_distribution_spec(1);
    x2=T_distribution_spec(2)+T_distribution_spec(1);
    x3=xlen-x2;
    x4=xlen-x1;

    slope=(T(2)-T(1))/xlen;

    coef2=[slope T(1)];
    
    T1=coef2(1)*x2+coef2(2);
    T2=coef2(1)*x3+coef2(2);
    
    M1=[x1^5 x1^4 x1^3 x1^2 x1 1
        x2^5 x2^4 x2^3 x2^2 x2 1
        5*x1^4 4*x1^3 3*x1^2 2*x1 1 0
        5*x2^4 4*x2^3 3*x2^2 2*x2 1 0
        20*x1^3 12*x1^2 6*x1 2 0 0
        20*x2^3 12*x2^2 6*x2 2 0 0];
    b1=[T(1);T1;0;slope;0;0];
    coef1=M1\b1;
    
    M2=[x3^5 x3^4 x3^3 x3^2 x3 1
        x4^5 x4^4 x4^3 x4^2 x4 1
        5*x3^4 4*x3^3 3*x3^2 2*x3 1 0
        5*x4^4 4*x4^3 3*x4^2 2*x4 1 0
        20*x3^3 12*x3^2 6*x3 2 0 0
        20*x4^3 12*x4^2 6*x4 2 0 0];
    b2=[T2;T(2);slope;0;0;0];
    coef3=M2\b2;
    
    
    for e=1:1:ne
        
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
%         yPos0=(iny-1)*dy;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            
            if xPos<=x1
                T_abs=T(1);
                dT_dx=0;
            elseif xPos<=x2
                T_abs=dot(coef1,[xPos^5 xPos^4 xPos^3 xPos^2 xPos 1]);
                dT_dx=dot(coef1,[5*xPos^4 4*xPos^3 3*xPos^2 2*xPos 1 0]);
            elseif xPos<=x3
                T_abs=dot(coef2,[xPos 1]);
                dT_dx=slope;
            elseif xPos<=x4
                T_abs=dot(coef3,[xPos^5 xPos^4 xPos^3 xPos^2 xPos 1]);
                dT_dx=dot(coef3,[5*xPos^4 4*xPos^3 3*xPos^2 2*xPos 1 0]);
            elseif xPos<=xlen
                T_abs=T(2);
                dT_dx=0;
            end
            
            
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
                temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
    
elseif distribution_type==8
    
    % 5th degree polynomial
    
    if length(T_distribution_spec)~=1 || T_distribution_spec>=xlen/2
        error('T_distribution_spec is the width of flat temperature.')
    end
    
    x1=T_distribution_spec;
    x2=xlen-x1;
    
    M1=[x1^5 x1^4 x1^3 x1^2 x1 1
        x2^5 x2^4 x2^3 x2^2 x2 1
        5*x1^4 4*x1^3 3*x1^2 2*x1 1 0
        5*x2^4 4*x2^3 3*x2^2 2*x2 1 0
        20*x1^3 12*x1^2 6*x1 2 0 0
        20*x2^3 12*x2^2 6*x2 2 0 0];
    b1=[T(1);T(2);0;0;0;0];
    coef=M1\b1;
    
%     coef
    
    for e=1:1:ne
        
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
%         yPos0=(iny-1)*dy;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            
            if xPos<=x1
                T_abs=T(1);
                dT_dx=0;
            elseif xPos<=x2
                T_abs=dot(coef,[xPos^5 xPos^4 xPos^3 xPos^2 xPos 1]);
                dT_dx=dot(coef,[5*xPos^4 4*xPos^3 3*xPos^2 2*xPos 1 0]);
            elseif xPos<=xlen
                T_abs=T(2);
                dT_dx=0;
            end
            
            
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
                temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
    
elseif distribution_type==9
    
    % Seventh degree polynomial
    
    if length(T_distribution_spec)~=1 || T_distribution_spec>=xlen/2
        error('T_distribution_spec is the width of flat temperature.')
    end
    
    x1=T_distribution_spec;
    x2=xlen-x1;
    
    M1=[x1^7 x1^6 x1^5 x1^4 x1^3 x1^2 x1 1
        x2^7 x2^6 x2^5 x2^4 x2^3 x2^2 x2 1
        7*x1^6 6*x1^5 5*x1^4 4*x1^3 3*x1^2 2*x1 1 0
        7*x2^6 6*x2^5 5*x2^4 4*x2^3 3*x2^2 2*x2 1 0
        42*x1^5 30*x1^4 20*x1^3 12*x1^2 6*x1 2 0 0
        42*x2^5 30*x2^4 20*x2^3 12*x2^2 6*x2 2 0 0
        210*x1^4 120*x1^3 60*x1^2 24*x1 6 0 0 0
        210*x2^4 120*x2^3 60*x2^2 24*x2 6 0 0 0];
    b1=[T(1);T(2);0;0;0;0;0;0];
    coef=M1\b1;
    
    
    for e=1:1:ne
        
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            
            if xPos<=x1
                T_abs=T(1);
                dT_dx=0;
            elseif xPos<=x2
                T_abs=dot(coef,[xPos^7 xPos^6 xPos^5 xPos^4 xPos^3 xPos^2 xPos 1]);
                dT_dx=dot(coef,[7*xPos^6 6*xPos^5 5*xPos^4 4*xPos^3 3*xPos^2 2*xPos 1 0]);
            elseif xPos<=xlen
                T_abs=T(2);
                dT_dx=0;
            end
            
            
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
                temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
	
elseif distribution_type==10
    
    % Polynomial with odd degree specified
    
    if length(T_distribution_spec)~=2 || T_distribution_spec(2)>=xlen/2
        error('T_distribution_spec=[degree, width of flat area]')
    end
    
    x1=T_distribution_spec(2);
    x2=xlen-x1;
    
    %-----------
%     warning('testing this line')
%     x1=0;
    %------------
    
    degree_plus=T_distribution_spec(1)+1;
    
    mx=mx_polynomial_order_constraint(T_distribution_spec(1),x1,x2);
    b=zeros(degree_plus,1);
    b(1)=T(1);
    b(2)=T(2);
    coef=mx\b;
    
    multiplier1=zeros(1,degree_plus);
    multiplier2=zeros(1,degree_plus);
    
    for e=1:1:ne
        
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            
            if xPos<x1
                T_abs=T(1);
                dT_dx=0;
            elseif xPos<=x2
                
                for i=1:1:degree_plus
                    multiplier1(i)=xPos^(degree_plus-i);
                end
                for i=1:1:T_distribution_spec(1)
                    
                    multiplier2(i)=(degree_plus-i)*xPos^(T_distribution_spec(1)-i);
                    
                    
                end
                
                T_abs=dot(coef,multiplier1);
                dT_dx=dot(coef,multiplier2);
            elseif xPos<=xlen
                T_abs=T(2);
                dT_dx=0;
            end
            
            
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
                temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
    
elseif distribution_type==11
    
    if length(T_distribution_spec)~=2 || T_distribution_spec(2)>=0.5
        error('T_distribution_spec=[degree, width of flat area]')
    end
    
    if length(n2)==2
        error('This distribution has not been designed for varying n2')
    end
    
    if chi_type~=1
        error('This is only desgined for chi_type==1')
    end
    
    error('Has not coded the temperature gradient yet')
    
    c2=1-ci_ave;
    
    chi1=chi_a+chi_b/T(1);
    k1=sqrt(-0.5*diff*T(1)*(1/(ci_ave*n1)+1/(c2*n2)-2*chi1/n1));
    
    chi2=chi_a+chi_b/T(2);
    k2=sqrt(-0.5*diff*T(2)*(1/(ci_ave*n1)+1/(c2*n2)-2*chi2/n1));
    
    x1=T_distribution_spec(2);
    x2=xlen-x1;
    degree_plus=T_distribution_spec(1)+1;
    mx=mx_polynomial_order_constraint(T_distribution_spec(1),x1,x2);
    b=zeros(degree_plus,1);
    b(1)=k1;
    b(2)=k2;
    coef=mx\b;
    
    multiplier1=zeros(1,degree_plus);
    multiplier2=zeros(1,degree_plus);
    
    a=-2/diff;
    b=2*chi_b/n1;
    c=1/(ci_ave*n1)+1/(c2*n2)-2*chi_a/n1;
    
    
    for e=1:1:ne
        
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            
            if xPos<x1
                T_abs=T(1);
%                 dT_dx=0;
%                 T_abs=k1;

%                 T_abs=(c1*k1^2+c2)/c3;

            elseif xPos<=x2
                
                for i=1:1:degree_plus
                    multiplier1(i)=xPos^(degree_plus-i);
                end
                for i=1:1:T_distribution_spec(1)
                    
                    multiplier2(i)=(degree_plus-i)*xPos^(T_distribution_spec(1)-i);
                end
                
                k=dot(coef,multiplier1);
                
                T_abs=(a*k^2+b)/c;
                
%                 T_abs=dot(coef,multiplier1);
%                 dT_dx=dot(coef,multiplier2);

%                 T_abs=k;
                
                
            elseif xPos<=xlen
                T_abs=T(2);
%                 dT_dx=0;

%                 T_abs=k2;

            end
            
            
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
%                 temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
    
elseif distribution_type==12
    
    error('This is not ready')
    
    if length(T_distribution_spec)~=2 || T_distribution_spec(2)>=0.5
        error('T_distribution_spec=[degree, width of flat area]')
    end
    
    if length(n2)==2
        error('This distribution has not been designed for varying n2')
    end
    
    if chi_type~=1
        error('This is only desgined for chi_type==1')
    end
    
    
    c2=1-ci_ave;
    
    chi1=chi_a+chi_b/T(1);
    k1=sqrt(-0.5*diff*T(1)*(1/(ci_ave*n1)+1/(c2*n2)-2*chi1/n1));
    R1=-k1^2*(diff*T(1)*(1/(ci_ave*n1)+1/(c2*n2)-2*chi1/n1)+k1^2);
    
    chi2=chi_a+chi_b/T(2);
    k2=sqrt(-0.5*diff*T(2)*(1/(ci_ave*n1)+1/(c2*n2)-2*chi2/n1));
    R2=-k2^2*(diff*T(2)*(1/(ci_ave*n1)+1/(c2*n2)-2*chi2/n1)+k2^2);
    
    x1=T_distribution_spec(2);
    x2=xlen-x1;
    degree_plus=T_distribution_spec(1)+1;
    mx=mx_polynomial_order_constraint(T_distribution_spec(1),x1,x2);
    b=zeros(degree_plus,1);
    b(1)=k1;
    b(2)=k2;
    coef=mx\b;
    
    multiplier1=zeros(1,degree_plus);
    multiplier2=zeros(1,degree_plus);
    
    a=-2/diff;
    b=2*chi_b/n1;
    c=1/(ci_ave*n1)+1/(c2*n2)-2*chi_a/n1;
    
    
%     R1
%     R2
% 
%     return
    
    
    for e=1:1:ne
        
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            
            if xPos<x1
                T_abs=T(1);
%                 dT_dx=0;
%                 T_abs=k1;

%                 T_abs=(c1*k1^2+c2)/c3;
                T_abs=R1;

            elseif xPos<=x2
                
                for i=1:1:degree_plus
                    multiplier1(i)=xPos^(degree_plus-i);
                end
                for i=1:1:T_distribution_spec(1)
                    
                    multiplier2(i)=(degree_plus-i)*xPos^(T_distribution_spec(1)-i);
                end
                
                k=dot(coef,multiplier1);
                
                T_abs=(a*k^2+b)/c;
                
%                 T_abs=dot(coef,multiplier1);
%                 dT_dx=dot(coef,multiplier2);

                chi=chi_a+chi_b/T_abs;

                T_abs=-k^2*(diff*T_abs*(1/(ci_ave*n1)+1/(c2*n2)-2*chi/n1)+k^2);
                
                
            elseif xPos<=xlen
                T_abs=T(2);
%                 dT_dx=0;

%                 T_abs=k2;
                T_abs=R2;

            end
            
            
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
%                 temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
    
elseif distribution_type==13
    
    %Sine wave across the domain
    
    if length(T_distribution_spec)~=1
        error('T_distribution_spec=number of periods within the domain')
    end
    
    mid_height=mean(T);
    amp=T(2)-mid_height;
    w=2*pi*T_distribution_spec/xlen;
    sin_shift=xlen/(4*T_distribution_spec);
    chunk=amp*w;
    
    for e=1:1:ne
        [inx,~]=inxiny_elemental(e,ney);
        xPos0=(inx-1)*dx;
        for ix=1:1:3
            xPos=xPos0+gp(ix);
            T_abs=amp*sin(w*(xPos-sin_shift))+mid_height;
            dT_dx=chunk*cos(w*(xPos-sin_shift));
            for iy=1:1:3
                temperature(e,ix,iy,1)=T_abs;
                temperature(e,ix,iy,2)=dT_dx;
            end
        end
    end
    
else
    error('This type of temperature distribution is not available')
end

end