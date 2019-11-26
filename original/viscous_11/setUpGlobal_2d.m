function [ne,nx,ny,n,nfour,ll,li,lu,bl,bi,bu,tl,ti,tu,rl,ri,ru,...
    fc,frequency_domain,x_coord,y_coord,T,...
    nnz_,irow,icol,weights,ci_ave,dxdy,...
    Two_chi_n1,coef_T,cp,bypass,dto,...
    nCoreTasks,index_array,...
    grad_T,...
    wTerms,iZeros,iOnes,k_target]=...
    setUpGlobal_2d(nex,ney,fc,alpha,xlen,ylen,...
    ci,T,entropy,diff,n1,n2,...
    communication_reliever,nworkers,...
    dt,distribution_type,T_distribution_spec,...
    export_figure,k_specific,...
    time_dependent_temperature,temperature_spec)

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

if length(T)==2 || time_dependent_temperature==1
    grad_T=1;
%     %T=linspace(T(1),T(2),nx);
%     if T(1)~=T(2)
%         grad_T=1;
%     elseif T(1)==T(2)
%         %-------------TEMPORARY----------------!!!!!!!!!
% %         T=T(1);
% %         grad_T=0;
% 
%         grad_T=1;
%         %----------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     end
elseif length(T)==1
%     T=ones(1,nx)*T;
    grad_T=0;
end

nnz_=nnzSJ(nx,ny);
[irow,icol]=setUpSparse(nnz_,ny,n);

weights=generateWeights();

dx=x_coord(2)-x_coord(1);
dy=y_coord(2)-y_coord(1);
weights=weight_adjuster(weights,dx,dy);

if time_dependent_temperature==0

    if grad_T==0
        coef_T=diff*T;
        Two_chi_n1=get_Two_chi_n1_with_order(T,entropy,n1,0);
    elseif grad_T==1
    %     [Two_chi_n1,coef_T]=chi_T_gradient_2d(T,dx,entropy,n1);


        coef_T=T_characterization(distribution_type,ne,T,...
            T_distribution_spec,nex,ney,xlen,ylen);
        Two_chi_n1=get_two_chi_n1(ne,coef_T,entropy,n1);
    end
    
elseif time_dependent_temperature==1
%     coef_T=NaN;
%     Two_chi_n1=NaN;
    coef_T=T_characterization(distribution_type,ne,[T(1) T(1)],...
        T_distribution_spec,nex,ney,xlen,ylen);
    Two_chi_n1=get_two_chi_n1(ne,coef_T,entropy,n1);
end


dxdy=dx*dy;

% if T(1)==T(nx)
%     T=T(1);
% else
%     T=[T(1) T(nx)];
% end

% if grad_T==1
% 
%     T=[T(1) T(nx)];
% end

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

% if grad_T==1
%     term_persistent=zeros(ne,3,3);
%     for e=1:1:ne
%         for ix=1:1:3
%             for iy=1:1:3
%                 term_persistent(e,ix,iy)=...
%                     get_Two_chi_n1_with_order(...
%                     coef_T(e,ix,iy,1),entropy,n1,1)...
%                     *coef_T(e,ix,iy,1)*2/n1;
%             end
%         end
%     end
% %     term_persistent=zeros(nex,3);
% %     dT_dx=(T(2)-T(1))/xlen;
% %     for i=1:1:nex
% %         term_persistent(i,:)=get_Two_chi_n1_with_order(coef_T(i,:),entropy,n1,1)*dT_dx;
% %     end
% %     %two_slope=2*dT_dx;
% elseif grad_T==0
%     term_persistent=NaN;%two_slope=NaN;
% end

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

end