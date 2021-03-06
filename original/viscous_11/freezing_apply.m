function [sf_new,sj_new]=freezing_apply(...
    type,spec,sf,sj,freezing_zeros,freezing_ones,...
    time,dt,xlen,dx,nx,ny,co,irow)

sf_new=sf;
sj_new=sj;

time_new=time+dt;

if type==1
    
    error('Need to fix this')
    
    if time_new<spec(1)
        return
    else
        
        freezing_duration=time_new-spec(1);
        
        if freezing_duration>=spec(2)
            
            freezing_duration
            warning('Im just stopping it with error for now')
            error('All have been frozen')
        else
            
            freezing_ratio=freezing_duration/spec(2);
            nBlocks=floor(freezing_ratio*xlen/dx)+1;
            
            freezing_ratio
            
            %-----------------------------
            [~,n]=size(freezing_zeros);
            for sub_block=1:1:n
                if freezing_zeros(1,sub_block)==0
                    break
                else
                    sj_new(freezing_zeros(1,sub_block))=0;
                end
            end
            for block=2:1:nBlocks
                for sub_block=1:1:n
                    sj_new(freezing_zeros(block,sub_block))=0;
                end
            end
            for isf=1:1:ny*4*nBlocks
                sf_new(isf)=0;
            end
            [~,n]=size(freezing_ones);
            for block=1:1:nBlocks
                for sub_block=1:1:n
                    sj_new(freezing_ones(block,sub_block))=1;
                end
            end
            
            %--------------------------------------------------------
            
        end
    end
    
    
elseif type==2
    
    error('This type is not ready')
    
%     acceleration=2*(xlen-spec(1)*spec(2))/spec(2)^2;
%     
%     acceleration
%     
%     time_new
% %     time_new=spec(2);
% %     time_new
%     
%     
%     distance=spec(1)*time_new+0.5*acceleration*time_new^2;
%     
%     
%     distance
%     
%     nBlocks=floor(distance/dx)+1;
%     
%     
%     nBlocks
%     
%     x=linspace(0,spec(2),100);
%     y=x;
%     for ijk=1:1:100
%         y(ijk)=spec(1)*x(ijk)+0.5*acceleration*x(ijk)^2;
%     end
%     plot(x,y)
%     
%     %------------------------
%     
% %     duration=10^-7;
%     
%     degree=3;
% %     x1=0;
% %     x2=xlen;
%     
%     t1=0;
%     t2=10^-7;
%     
%     vo=xlen/10^-10;
%     vf=xlen/10^-5;
%     
%     mx=mx_polynomial_order_constraint(degree,t1,t2);
%     
%     constraints=zeros(degree+1,1);
%     constraints(2)=xlen;
%     constraints(3)=vo;
%     constraints(4)=vf;
%     
%     size(mx)
%     size(constraints)
%     
%     coefs=mx\constraints;
%     
%     coefs=coefs'
%     
%     coefs
%     
%     
%     x=linspace(t1,t2,100);
%     y=x;
%     for ijk=1:1:100
%         y(ijk)=dot(coefs,[x(ijk)^3 x(ijk)^2 x(ijk) x(ijk)]);
%     end
%     plot(x,y)
    
    
%     spec=[]

    
    vo=(xlen-0.5*spec(2)*spec(1)^2)/spec(1);
    
    x=linspace(0,spec(1),100);
    y=x;
    for ijk=1:1:100
        
        0.5*spec(2)*x(ijk)^2
        
        y(ijk)=vo*x(ijk)+0.5*spec(2)*x(ijk)^2;
    end
    plot(x,y)
    
    
    error('just stop here')
    
    
    %-----------------------------
    [~,n]=size(freezing_zeros);
    for sub_block=1:1:n
        if freezing_zeros(1,sub_block)==0
            break
        else
            sj_new(freezing_zeros(1,sub_block))=0;
        end
    end
    for block=2:1:nBlocks
        for sub_block=1:1:n
            sj_new(freezing_zeros(block,sub_block))=0;
        end
    end
    for isf=1:1:ny*4*nBlocks
        sf_new(isf)=0;
    end
    [~,n]=size(freezing_ones);
    for block=1:1:nBlocks
        for sub_block=1:1:n
            sj_new(freezing_ones(block,sub_block))=1;
        end
    end

    %--------------------------------------------------------
    
elseif type==3
    
    
    c_nodal=extract_nodal_weights_2D(co,nx,ny);
    
    x_coord=linspace(0,xlen,nx);
    trend=zeros(1,nx);
    
    
    
    for xth=1:1:nx
        summation=0;
        mean_val=mean(c_nodal(:,xth));
        for yth=1:1:ny
            summation=summation+(c_nodal(yth,xth)-mean_val)^2;
        end
        trend(xth)=summation/ny;
    end
    
    plot(x_coord,trend)
    drawnow
    
    
    
%     nBlocks=0;
    
    frozen_blocks=zeros(1,nx);
    
    for xth=1:1:nx
        if trend(xth)>=spec
            frozen_blocks(xth)=1;
        end
    end
    
    nBlocks=sum(frozen_blocks);
    
    
    size(freezing_ones)
    nnz(freezing_ones)
    
    
    
    if nBlocks==0
        return
    elseif nBlocks==nx
        sf_new=NaN;
        sj_new=NaN;
        return
    end
    
    %-----------------------------
    
    [sf_new,sj_new]=applying_freezing(...
        sf_new,sj_new,freezing_zeros,freezing_ones,...
        frozen_blocks,irow,nx);
    
    frozen_blocks'
    
%     error('yaminoma')
    %--------------------------------------------------------
    
%     size(c_nodal)
%     
%     
    

    
else
    error('Specified type does not exist')
end


end

function [sf_new,sj_new]=applying_freezing(...
    sf_new,sj_new,freezing_zeros,freezing_ones,...
    frozen_blocks,irow,nx)

[~,n]=size(freezing_zeros);

if frozen_blocks(1)==1
    for sub_block=1:1:n
        if freezing_zeros(1,sub_block)==0
            break
        else
            sj_new(freezing_zeros(1,sub_block))=0;
        end
    end
end


for block=2:1:length(frozen_blocks)-1
    if frozen_blocks(block)==1
        for sub_block=1:1:n
            sj_new(freezing_zeros(block,sub_block))=0;
        end
    end
end


if frozen_blocks(nx)==1
    for sub_block=1:1:n
        if freezing_zeros(1,sub_block)==0
            break
        else
            sj_new(freezing_zeros(1,sub_block))=0;
        end
    end
end

%---------

[m,n]=size(freezing_ones);

for block=1:1:m
    if frozen_blocks(block)==1
        for sub_block=1:1:n
            sj_new(freezing_ones(block,sub_block))=1;
            sf_new(irow(freezing_ones(block,sub_block)))=0;
        end
    end
end

%---------


% return
% 
% %-----------------------------
% 
% for isf=1:1:ny*4*nBlocks
%     sf_new(isf)=0;
% end
% [~,n]=size(freezing_ones);
% for block=1:1:nBlocks
%     for sub_block=1:1:n
%         sj_new(freezing_ones(block,sub_block))=1;
%         
%     end
% end



end