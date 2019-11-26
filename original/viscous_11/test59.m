function test59

clear
clc


eX=[1 5 2 7 0 4 1 3];
eY=[-1 0 1 2 5 4 5 6];
eZ=[0 -1 3 4 -2 0 4 6];

% hold on
% for i=1:1:8
%     plot3(eX(i),eY(i),eZ(i),'o')
% end
% hold off
% grid on

connectivity=[...
    1 2
    3 4
    1 3
    2 4
    5 6
    7 8
    5 7
    6 8
    1 5
    2 6
    3 7
    4 8];
% hold on
% for i=1:1:12
%     plot3([eX(connectivity(i,1)),eX(connectivity(i,2))],[eY(connectivity(i,1)),eY(connectivity(i,2))],[eZ(connectivity(i,1)),eZ(connectivity(i,2))])
% end
% hold off

c=rand(64,1)*10;
what_to_obtain=1;

alpha=0.2;
beta=0.3;
gamma=0.4;

x=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 0 0]);
y=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 0 0]);
z=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,[0 0 0]);

% hold on
% plot3(x,y,z,'rx')
% hold off

ws=[0 0 0];
ts=[0 0 0];

orientation_list=...
    [0 0 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1
    1 0 1
    0 1 1
    1 1 1];

type_list=...
    [0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    1 0 1
    0 1 1
    1 1 1];
    

%--------------------
invJ=get_invJ(alpha,beta,gamma,eX,eY,eZ);

J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ);
inv(J1);

f_local_derivatives=get_f_dervs_3d(alpha,beta,gamma,ws,ts);
test_original=get_fg_dervs_3d(alpha,beta,gamma,ws,ts,eX,eY,eZ,what_to_obtain);
test1=0;

ic=0;
for iorientation=1:1:8
    for itype=1:1:8
        ic=ic+1;
        test1=test1+...
            c(ic)*get_fg_dervs_3d(alpha,beta,gamma,orientation_list(iorientation,:),type_list(itype,:),eX,eY,eZ,what_to_obtain);
    end
end
test1;
%--------------------
test2=0;
ic=0;
for iorientation=1:1:8
    for itype=1:1:8
        ic=ic+1;
        f_local_derivatives=get_f_dervs_3d(alpha,beta,gamma,orientation_list(iorientation,:),type_list(itype,:));
        
        for ijk=1:1:9
            f_local_derivatives(ijk)=f_local_derivatives(ijk)*c(ic);
        end
        
        take=0;
        for ijk=1:1:9
            take=take+invJ(what_to_obtain,ijk)*f_local_derivatives(ijk);
        end
        test2=test2+take;
        
        
    end
end


test2;

test1-test2;

%------------------------------------------------

alpha=1;
beta=0;
gamma=1;

eX=[0 1 0 1 0 2 0 1];
eY=[0 0 1 1 0 0 1 1];
eZ=[0 0 0 0 1 1 1 1];

map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 0 0]);

J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ);

%-----------------------------------------------------------


% dir1=[1 1 -2];
% dir2=[1 2 -1];
% dir3=[2 -2 1];
% dir4=[1 3 3];
% 
% dir1=make_unit(dir1)

dirs=[...
    1 1 -2
	1 2 -1
	2 -2 1
	1 3 3];

len=0.8;

for i=1:1:4
    take=make_unit(dirs(i,:));
    eX(i)=take(1);
    eY(i)=take(2);
    eZ(i)=take(3);
    
    
    eX(4+i)=take(1)*len;
    eY(4+i)=take(2)*len;
    eZ(4+i)=take(3)*len;

%     eX(4+i)=0;
%     eY(4+i)=0;
%     eZ(4+i)=0;
    
end

hold on
for i=1:1:12
    plot3([eX(connectivity(i,1)),eX(connectivity(i,2))],[eY(connectivity(i,1)),eY(connectivity(i,2))],[eZ(connectivity(i,1)),eZ(connectivity(i,2))])
end
hold off
grid on

hold on

take=4;
plot3(eX(take),eY(take),eZ(take),'rx')
hold off

alpha=0.3;
beta=0.5;
gamma=0;

J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ);

take=J1(3,:);
take101=take;
sqrt(sum(take.^2));

take=make_unit(take);

xyz1=map_to_global_3_compts_3d(alpha,beta,gamma,eX,eY,eZ,[0 0 0]);
xyz2=map_to_global_3_compts_3d(alpha,beta,1,eX,eY,eZ,[0 0 0]);

take2=make_unit(xyz1-xyz2);

take+take2;

invJ=inv(J1);

what_to_obtain=[1 2 3];
nObtaining=length(what_to_obtain);

test1=zeros(1,nObtaining);
ic=0;
for iorientation=1:1:8
    for itype=1:1:8
        take=zeros(1,nObtaining);
        ic=ic+1;
        f_local_derivatives=get_f_dervs_3d(alpha,beta,gamma,orientation_list(iorientation,:),type_list(itype,:));
        
        for ijk=1:1:9
%             [iorientation itype ic]
            f_local_derivatives(ijk)=f_local_derivatives(ijk)*c(ic);
        end
        if type_list(itype,3)==1
            f_local_derivatives(3)=0;
        end
        
        for ij=1:1:nObtaining
            for ijk=1:1:3
%                 ij
%                 ijk
%                 what_to_obtain(ij)
%                 size
%                 invJ(what_to_obtain(ij),ijk)
                take(ij)=take(ij)+invJ(what_to_obtain(ij),ijk)*f_local_derivatives(ijk);
            end
        end
        
        test1=test1+take;
        
    end
end

test1;
dot(test1,take101);

%----------------------------------------------------------

what_to_obtain=[1 2 3];
nObtaining=length(what_to_obtain);

test1=zeros(1,nObtaining);
ic=0;
test111=0;
invJ=get_invJ(alpha,beta,gamma,eX,eY,eZ);
for iorientation=1:1:8
    for itype=1:1:8
        take=zeros(1,nObtaining);
        ic=ic+1;
        f_local_derivatives=get_f_dervs_3d(alpha,beta,gamma,orientation_list(iorientation,:),type_list(itype,:));
        
%         for ijk=1:1:9
% %             [iorientation itype ic]
%             f_local_derivatives(ijk)=f_local_derivatives(ijk)*c(ic);
%         end
%         if type_list(itype,3)==1
%             f_local_derivatives(3)=0;
%         end
%         
%         [sum(f_local_derivatives(1:1:3)) iorientation itype]

        
        
        for ij=1:1:nObtaining
            for ijk=1:1:9
%                 ij
%                 ijk
%                 what_to_obtain(ij)
%                 size
%                 invJ(what_to_obtain(ij),ijk)
                take(ij)=take(ij)+invJ(what_to_obtain(ij),ijk)*f_local_derivatives(ijk);
            end
        end
        
%         sum(f_local_derivatives(1:1:3));
%         [sum(f_local_derivatives(1:1:3)) iorientation itype gamma];
        
%         if type_list(itype,3)==1 && ic<=64
%             c(ic)=0;
%             test111=test111+1;
%             [test111 iorientation itype sum(f_local_derivatives(1:1:3))];
%         end
        
        test1=test1+take*c(ic);
        
    end
end

test1;
dot(test1,take101);

[alpha beta gamma];
f_local_derivatives=get_f_dervs_3d(alpha,beta,gamma,[1 0 1],[0 0 0]);
f_local_derivatives(1:1:3);

%--------------------------------------------------

% figure(2)

eX=[2 5 2 7 0 4 1 3];
eY=[-1 0 1 2 5 4 5 6];
eZ=[0 -1 3 4 -2 0 4 6];


% eX=[0 1 0 1 0 2 0 1];
% eY=[0 0 1 1 0 0 1 1];
% eZ=[0 0 0 0 1 1 1 1];

% eX=[0 2 0 1 0 2 0 2];
% eY=[0 0 1 1 0 0 1 1];
% eZ=[0 0 0 0 1 1 1 1];

c;
alpha=0.2;
beta=0.3;
gamma=0.6;
invJ=get_invJ(alpha,beta,gamma,eX,eY,eZ);
for i=1
    plot3([eX(connectivity(i,1)),eX(connectivity(i,2))],[eY(connectivity(i,1)),eY(connectivity(i,2))],[eZ(connectivity(i,1)),eZ(connectivity(i,2))])
end
hold on
for i=2:1:12
    plot3([eX(connectivity(i,1)),eX(connectivity(i,2))],[eY(connectivity(i,1)),eY(connectivity(i,2))],[eZ(connectivity(i,1)),eZ(connectivity(i,2))])
end
hold off
grid on

c=rand(1,64);
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2]
w=[5/18 4/9 5/18];
[cs,weights]=compute_c_specific(c,1:1:7,0,0,0,gps,w,eX,eY,eZ,orientation_list);

cs'
c(1:1:8)

orientation_list;

disp('hey---------------------')


weights(7,2,2);

for orientation=1:1:8
    for type=1:1:8
        for order=3
            if weights(orientation,type,order)~=0
                [orientation,type,order,weights(orientation,type,order)]
            end
        end
    end
end

disp('done---------------------')

disp('confirmation')
c(9:1:16)
cs(1:1:4)
c(1:1:4);

disp('2d test--------------------')

eX=[1 5 0 7];
eY=[-2 -3 1 2];

% plot(eX,eY,'*')
% hold on
% node=1;
% plot(eX(node),eY(node),'rx')
% hold off
% grid on
% xlabel('x');ylabel('y')

mx=compute_c_specific_2d(...
    1,1,eX,eY);

take=temp_dim_reducer(mx,3);


disp('MONDAY---------------------------------------')

eX=[2 5 2 7 0 4 1 3];
eY=[-1 0 1 2 5 4 5 6];
eZ=[0 -1 3 4 -2 0 4 6];

coeffs=get_Hermitian_pol_coeffs_3d();

alpha=1;
beta=0;
gamma=1;
order=1;

c=rand(1,32);

mx=1:1:32;

mx=[mx' c']

% c=ones(1,32);

take=compute_specific_c_3d(c,alpha,beta,gamma,eX,eY,eZ,coeffs);

take

phis=get_phis();
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

weights=compute_elemental_weights_serendipity(gps,phis,eX,eY,eZ);

for ix=1
    for iy=2
        for iz=3
            
        end
    end
end

%-------------------------------------------------------

size(weights)

% c=ones(1,32);

c'

conc_=compute_conc_serendipity_(c,weights);

ix=1;
iy=2;
iz=3;
order=2;

% gps=zeros(1,3);

take=compute_specific_c_3d(c,gps(ix),gps(iy),gps(iz),eX,eY,eZ,coeffs);

conc_(:,ix,iy,iz)
take

%----------------------------------------------------------

clear
clc

disp('Thursday--------------------------')

eX=[2 5 2 7 0 4 1 3];
eY=[-1 0 1 2 5 4 5 6];
eZ=[0 -1 3 4 -2 0 4 6];

coeffs=get_Hermitian_pol_coeffs_3d();
c=rand(1,32);
gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
phis=get_phis();
weights=compute_elemental_weights_serendipity(gps,phis,eX,eY,eZ);

conc_=compute_conc_serendipity_(c,weights);
summation=0;

for ix=1:1:3
    for iy=1:1:3
        for iz=1:1:3
            take=compute_specific_c_3d(c,gps(ix),gps(iy),gps(iz),eX,eY,eZ,coeffs);
            for order=1:1:7
                
                diff=abs(take(order)-conc_(order,ix,iy,iz));
                summation=summation+diff;
            end
        end
    end
end

summation



end

%======================================
%======================================
%======================================
%======================================
%======================================
%======================================
%Accessory functions=====================================
%======================================
%======================================
%======================================
%======================================
%======================================
%======================================

function sol=temp_dim_reducer(mx,col3rd)
sol=zeros(4,4);
for dim1=1:1:4
    for dim2=1:1:4
        sol(dim1,dim2)=mx(dim1,dim2,col3rd);
    end
end
end

function sol=make_unit(vector)
len=sqrt(sum(vector.^2));
sol=vector/len;
end

%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

function J1=get_Jacobian_part_1_3d(alpha,beta,gamma,eX,eY,eZ)
J1=zeros(3,3);
order_list=...
    [1 0 0
    0 1 0
    0 0 1];
for i=1:1:3
    J1(i,1)=map_to_global_1_compt_3d(alpha,beta,gamma,eX,order_list(i,:));
end
for i=1:1:3
    J1(i,2)=map_to_global_1_compt_3d(alpha,beta,gamma,eY,order_list(i,:));
end
for i=1:1:3
    J1(i,3)=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,order_list(i,:));
end
end