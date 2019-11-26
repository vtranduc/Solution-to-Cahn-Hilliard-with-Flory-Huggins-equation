function cp=predict_c_with_acceleration_3d(c,co,coo,dt,dto,dtoo,...
    parallel_computing,n_gbfs)

cp=zeros(1,n_gbfs);

if parallel_computing==1
    parfor i=1:1:n_gbfs
        v2=(c(i)-co(i))/dto;
        cp(i)=...
            c(i)+v2*dt...
            +(v2-(co(i)-coo(i))/dtoo)/(dto+dtoo)*dt^2.0; 
    end
elseif parallel_computing==0
    v2=(c-co)/dto;
    cp=...
        c+v2*dt...
        +(v2-(co-coo)/dtoo)/(dto+dtoo)*dt^2.0; 
end

end