function [frequency_domain,magnitude]=fourier_analysis_3d(c_,ci_ave,alpha,beta,gamma,fc,nworkers)

magnitude=zeros(1,alpha);

[M,N,O]=size(c_);

c=c_-ci_ave*ones(M,N,O);

[allocation,extra,load_less,load_more]=wise_task_splitter_v2(alpha,nworkers);

temp=zeros(nworkers,load_more);

frequency_domain=linspace(0,fc,alpha);

twopi=2*pi;

angle_list_1=linspace(0,twopi-twopi/beta,beta);
angle_list_2=linspace(0,twopi-twopi/gamma,gamma);

nEvals=beta*gamma;

mx_angle=zeros(beta,gamma,3);
for ibeta=1:1:beta
    for igamma=1:1:gamma
        mx_angle(ibeta,igamma,1)=sin(angle_list_1(ibeta))*cos(angle_list_2(igamma));
        mx_angle(ibeta,igamma,2)=sin(angle_list_1(ibeta))*sin(angle_list_2(igamma));
        mx_angle(ibeta,igamma,3)=cos(angle_list_1(ibeta));
    end
end

parfor worker=1:1:nworkers
    if worker<=extra
        temp(worker,:)=fourier_analysis_3d_assist(allocation(worker,:),load_more,load_more,...
            c,frequency_domain,...
            mx_angle,nEvals,beta,gamma,M,N,O,twopi);
    elseif worker>extra
        temp(worker,:)=fourier_analysis_3d_assist(allocation(worker,:),load_less,load_more,...
            c,frequency_domain,...
            mx_angle,nEvals,beta,gamma,M,N,O,twopi);
    end
end

index=0;
for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    else
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        magnitude(index)=temp(worker,i);
    end
end

end

function sol=fourier_analysis_3d_assist(keys,load,load_more,...
    c,frequency_domain,...
    mx_angle,nEvals,beta,gamma,M,N,O,twopi)

sol=zeros(1,load_more);
for i=1:1:load
    k=frequency_domain(keys(i));
    summation=0;
    for ibeta=1:1:beta
        for igamma=1:1:gamma
            summation=summation...
                +get_mag_3d(c,...
                k*mx_angle(ibeta,igamma,1),...
                k*mx_angle(ibeta,igamma,2),...
                k*mx_angle(ibeta,igamma,3),...
                M,N,O,twopi);
        end
    end
    sol(i)=summation/nEvals;
end
end