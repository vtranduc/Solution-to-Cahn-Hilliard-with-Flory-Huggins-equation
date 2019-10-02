function [ne,nx,ny,n,nfour,ll,li,lu,bl,bi,bu,tl,ti,tu,rl,ri,ru,...
    fc,frequency_domain,x_coord,y_coord,T,...
    nnz_,irow,icol,weights,ci_ave,dxdy,...
    chi,coef_T,cp,bypass,dto,...
    ...
    ...
    wTerms,iZeros,iOnes,k_target,...
    coef_n2,wwTerms]=...
    setUpGlobal_2d(nex,ney,fc,alpha,xlen,ylen,...
    ci,T,entropy,diff,n1,n2,...
    communication_reliever,nworkers,...
    dt,distribution_type,T_distribution_spec,...
    export_figure,k_specific,...
    chi_a,chi_b,chi_type,...
    n2_distribution_type,n2_distribution_spec)

%======================
% RENAME coef_T IN ENTIRE PROGRAM!!!


ne=nex*ney;
nx=nex+1;ny=ney+1;
n=nx*ny;
nfour=4*n;
ll=2;li=4;lu=2+4*ney;
bl=3;bi=4*ny;bu=3+4*nex*ny;
tl=3+4*ney;ti=4*ny;tu=3+4*ney+4*nex*ny;
rl=2+4*nex*ny;ri=4;ru=2+4*ney+4*nex*ny;
if fc==0
    fc=1/(2*min([xlen/nex ylen/ney]));
end
frequency_domain=linspace(0,fc,alpha);

x_coord=linspace(0,xlen,nx);
y_coord=linspace(0,ylen,ny);

if length(ci)==1
    ci_ave=ci;
else
    ci_ave=(ci(1)+ci(2))/2;
end

% if length(T)==2 || time_dependent_temperature==1
%     grad_T=1;
% elseif length(T)==1
%     grad_T=0;
% end

nnz_=nnzSJ(nx,ny);
[irow,icol]=setUpSparse(nnz_,ny,n);

weights=generateWeights();

dx=x_coord(2)-x_coord(1);
dy=y_coord(2)-y_coord(1);
weights=weight_adjuster(weights,dx,dy);

%----Temperature---------
if length(T)==1
    coef_T=zeros(ne,3,3,3);
    for e=1:1:ne
        for ix=1:1:3
            for iy=1:1:3
                coef_T(e,ix,iy,1)=T;
            end
        end
    end
    T=[T T];
elseif length(T)==2
    coef_T=T_characterization(distribution_type,ne,T,...
    T_distribution_spec,nex,ney,xlen,ylen,...
    diff,n1,n2,ci_ave,chi_a,chi_b,chi_type);
end
chi=interaction_parameter_c_independent(coef_T,chi_a,chi_b,chi_type);



%---Polymerization-------------------

if length(n2)==1
    coef_n2=zeros(ne,3,3,3);
    for e=1:1:ne
        for ix=1:1:3
            for iy=1:1:3
                coef_n2(e,ix,iy,1)=n2;
            end
        end
    end
    
elseif length(n2)==2
    coef_n2=n2_characterization(n2_distribution_type,ne,n2,...
        n2_distribution_spec,nex,ney,xlen,ylen);
end

%------------------------------

dxdy=dx*dy;
cp=zeros(1,nfour);
bypass=0;
dto=dt;

if communication_reliever==1
    remainder=mod(nnz_,nworkers);
    if remainder==0
        njobs=nnz_/nworkers;
        nCoreTasks=njobs;
    else
        njobs=(nnz_-remainder)/nworkers;
        nCoreTasks=njobs+1;
    end
    index_array=zeros(nworkers,2);
    index_array(1,:)=[1 njobs];
    for job=2:1:nworkers
        if remainder~=0
            index_array(job-1,2)=index_array(job-1,2)+1;
            remainder=remainder-1;
        end
        index_array(job,1)=index_array(job-1,2)+1;
        index_array(job,2)=index_array(job-1,2)+njobs;
    end
elseif communication_reliever==0
    nCoreTasks=NaN;
    index_array=NaN;
end

wTerms=compute_wTerms_2d(weights);
[iZeros,iOnes]=get_bc_keyVals(nx,ny,nnz_,irow,icol);
if export_figure==1 && k_specific==0
    if length(T)>1 || length(ci)>1
        warning('k_target will be characteristic frequency of mean of temperature ends and mean of concentration ends, not the average over the volume!')
    end
    k_target=get_characteristic_frequency(diff,ci_ave,mean(T),n1,n2,entropy);
else
    k_target=NaN;
end

wwTerms=zeros(3,3,16,16,2);
for ilocal=1:1:16
    for jlocal=1:1:16
        for ix=1:1:3
            for iy=1:1:3
                wwTerms(ix,iy,ilocal,jlocal,1)=...
                    weights(ix,iy,ilocal,2)*weights(ix,iy,jlocal,2)...
                    +weights(ix,iy,ilocal,3)*weights(ix,iy,jlocal,3);
                wwTerms(ix,iy,ilocal,jlocal,2)=...
                    wTerms(ix,iy,ilocal)*wTerms(ix,iy,jlocal);
            end
        end
    end
end

end