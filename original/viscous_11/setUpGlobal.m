function [ne,nx,ny,n,nfour,ll,li,lu,bl,bi,bu,tl,ti,tu,rl,ri,ru,...
    fc,frequency_domain,isf1d,x_coord,y_coord,T,co,...
    nnz_,irow,icol,weights]=...
    setUpGlobal(nex,ney,fc,alpha,num_str_fac_1d,xlen,ylen,...
    T,ci,include_fluc,ci_fluc)
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

if num_str_fac_1d>nx
    error('num_str_fac_1d must be smaller than or equal to nx!')
end
isf1d=zeros(1,num_str_fac_1d);
spacer=nx/(num_str_fac_1d-1);
curr_pos=0;
for i=2:1:num_str_fac_1d-1
    curr_pos=curr_pos+spacer;
    isf1d(i)=round(curr_pos);
end
isf1d(1)=1;
isf1d(num_str_fac_1d)=nx;


% x=zeros(1,n); y=zeros(1,n);
% 
% for i=1:n
%     x(i)=(xlen/nex)*floor((i-1)/ny);
% end
% for i=1:ny
%     for j=1:nx
%         y(i+(j-1)*ny)=(ylen/ney)*(i-1);
%     end
% end

x_coord=linspace(0,xlen,nx);
y_coord=linspace(0,ylen,ny);

if length(T)==2
    T=linspace(T(1),T(2),nx);
elseif length(T)==1
    T=ones(1,nx)*T;
end
% chi=get_chi(T,entropy,T_theta);

co=generate_co_2D(ci,nx,ny,include_fluc,ci_fluc);
co=icc(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);

nnz_=nnzSJ(nx,ny);
[irow,icol]=setUpSparse(nnz_,ny,n);

weights=generateWeights();

end