function co=generate_co_2D(conc_,nx,ny,include_fluc,ci_fluc)

%conc_ must be either a scalar or a vector of 2 elements

[nrows,ncols]=size(conc_);
if (nrows==1 && ncols==2) || (nrows==2 && ncols==1)
    low_c=conc_(1);
    high_c=conc_(2);
    if include_fluc==1
        slope=(high_c-low_c)/(nx-1);
        co_nodal=zeros(ny,nx);
        for x_coord=1:1:nx-1
            c_specific=slope*(x_coord-1)+low_c;
            for y_coord=1:1:ny
                co_nodal(y_coord,x_coord)=c_specific+(ci_fluc*(2*rand(1,1)-1));
            end
            ci_diff=c_specific-mean(co_nodal(:,x_coord));
            for y_coord=1:1:ny
                co_nodal(y_coord,x_coord)=co_nodal(y_coord,x_coord)+ci_diff;
            end
        end
        for y_coord=1:1:ny
            co_nodal(y_coord,nx)=high_c+(ci_fluc*(2*rand(1,1)-1));
        end
        ci_diff=high_c-mean(co_nodal(:,nx));
        for y_coord=1:1:ny
            co_nodal(y_coord,nx)=co_nodal(y_coord,nx)+ci_diff;
        end
        co_nodal=reshape(co_nodal,[nx*ny,1]);
        co=zeros(1,4*nx*ny);
        for i=1:1:nx*ny
            co((i-1)*4+1)=co_nodal(i);
        end
    elseif include_fluc==0
        slope=(high_c-low_c)/(nx-1);
        co_nodal=zeros(ny,nx);
        for x_coord=1:1:nx-1
            c_specific=slope*(x_coord-1)+low_c;
            for y_coord=1:1:ny
                co_nodal(y_coord,x_coord)=c_specific;
            end
        end
        for y_coord=1:1:ny
            co_nodal(y_coord,nx)=high_c;
        end
        co_nodal=reshape(co_nodal,[nx*ny,1]);
        co=zeros(1,4*nx*ny);
        for i=1:1:nx*ny
            co((i-1)*4+1)=co_nodal(i);
        end
    else
        error('Include_fluc must be 0 or 1!')
    end
elseif nrows==1 && ncols==1
    if include_fluc==1
        co_nodal=zeros(1,nx*ny);
        for i=1:1:nx*ny
            co_nodal(i)=conc_+(ci_fluc*(2*rand(1,1)-1));
        end
        conc__mean_diff=conc_-mean(co_nodal);
        co=zeros(1,4*nx*ny);
        for i=1:1:nx*ny
            co(4*(i-1)+1)=co_nodal(i)+conc__mean_diff;
        end
    elseif include_fluc==0
        co=zeros(1,4*nx*ny);
        for i=1:1:nx*ny
            co(4*(i-1)+1)=conc_;
        end
    else
        error('Include_fluc must be 0 or 1!')
    end
else
    error('conc_ must be a scalar or vector of 2 elements!')
end

end