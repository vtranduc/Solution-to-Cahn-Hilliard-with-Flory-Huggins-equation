function cp=predict_c_3d(c,co,dt,dto,parallel_computing,n_gbfs)

cp=zeros(1,n_gbfs);

if parallel_computing==1
    parfor i=1:1:n_gbfs
        cp(i)=c(i)+dt*((c(i)-co(i))/dto);
    end
elseif parallel_computing==0
    cp=c+dt*(c-co)/dto;
end

end