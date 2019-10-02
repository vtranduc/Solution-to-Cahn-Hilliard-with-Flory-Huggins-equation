function co=generate_co_anisotropic_2d(ci,type,spec,n,nx,ny,nfour,ci_fluc,...
    xlen,ylen)

% Note that this function uses random normal distribution

co=zeros(1,nfour);

if length(ci)==1
    
    for i=1:4:nfour-3
        co(i)=(ci_fluc*(2*randn(1,1)-1))+ci;
    end
    
elseif length(ci)==2
    
    % Linear
    
    dx=xlen/(nx-1);   
    if type==1
        slope=(ci(2)-ci(1))/xlen;
        slope_dx=slope*dx;
        multiplier=0;
        index=-3;
        slope_dx_multiplier=slope_dx*multiplier;
        for i=1:1:n
            index=index+4;
            co(index)=ci(1)+slope_dx_multiplier+(ci_fluc*(2*randn(1,1)-1));
            co(index+1)=slope;
            if mod(i,ny)==0
                multiplier=multiplier+1;
                slope_dx_multiplier=slope_dx*multiplier;
            end
        end
    elseif type==2
        
        % third degree polynomial
        
        pos2=xlen-spec;
        
        M=[spec^3 spec^2 spec 1
            pos2^3 pos2^2 pos2 1
            3*spec^2 2*spec 1 0
            3*pos2^2 2*pos2 1 0];
        divisor=[ci(1);ci(2);0;0];
        coef=M\divisor;
        
        dx=xlen/(nx-1);
        distance=0;
        c_local=ci(1);
        slope_local=0;
        index=-3;
        for i=1:1:n
            index=index+4;
            co(index)=c_local+(ci_fluc*(2*randn(1,1)-1));
            co(index+1)=slope_local;
            if mod(i,ny)==0
                distance=distance+dx;
                
                if distance>spec
                    if distance<pos2
                        c_local=dot(coef,[distance^3 distance^2 distance 1]);
                        slope_local=dot(coef,[3*distance^2 2*distance 1 0]);
                    elseif distance>=pos2
                        c_local=ci(2);
                        slope_local=0;
                    end
                end
            end
            
        end
        
        %--------------------------Find out specific point
%         syms x
%         eqn=coef(1)*x^3+coef(2)*x^2+coef(3)*x+coef(4)==0.3054;
%         solve(eqn, x)
%         error('just stop here')
        %---------------------------
        
    elseif type==3
        
        % Sinusoidal. T2 at the center
        
        error('Still under construction')
        
        X=linspace(0,xlen,nx);
        Y=linspace(0,ylen,ny);
        
        for node=1:1:n
            
            yth=mod(node,ny);
            if yth==0
                yth=ny;
            end
            xth=(node-yth)/ny+1;
            index=(node-1)*4+1;
            
            
            co(index)=xth*yth;

            
        end
        
        
    else
        error('No such type is available')
    end
else
    error('The length of ci must be 1 or 2')
end


%-------------Adjust----------
ci_ave=mean(ci);
[co,~]=c_adjuster_2d(co,ci_ave,n,nfour);

%----------------------------


end