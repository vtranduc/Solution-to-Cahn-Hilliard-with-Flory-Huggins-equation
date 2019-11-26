function hey=get_inverse_Jacobian_3d(ne,nex,nx,ny,nexney,X,Y,Z)

%----- Real thing starts here

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

gps_plus=zeros(3,3);
for i=1:1:3
    gps_plus(i,:)=gps*gps(i);
end
% 
% gps_plus
% gps
% error('afdafd')
%--------------------------------

w=[5/18 4/9 5/18];

invJ_mx=zeros(ne,3,3,3,9,9);

weights_plus=generateWeights_plus_3d();

% f_dervs=get_f_dervs_plus_3d(gps_plus);

f_dervs=2;

for e=1 %PARFOR MUST BE FEASIBLE!!!
    [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z);
    invJ=e_invJ(gps,eX,eY,eZ);
    hey=map_weights_with_inverse_Jacobian_3d(invJ,weights_plus,f_dervs,gps,gps_plus,w,eX,eY,eZ);
end
%--------------------More advanced test

hey;

end

function sol=map_weights_with_inverse_Jacobian_3d(invJ,weights_plus,f_dervs,gps,gps_plus,w,eX,eY,eZ)
sol=zeros(8,8,7,3,3,3);


%   MUST REVIEW THE IMPORTANCE OF WEIGHT_PLUS!!!

% sol=get_fg_dervs_3d(alpha,beta,gamma,orientations,types,eX,eY,eZ)
test=zeros(8,3,3,3,9);
orientation_list=...
    [0 0 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1
    1 0 1
    0 1 1
    1 1 1];
types=...
    [0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    1 0 1
    0 1 1
    1 1 1];

% f_der=get_f_dervs_3d(alpha,beta,gamma,orientations,types)
for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
%                 test(iorientation,ix,iy,iz,:)=get_fg_dervs_3d(gps(ix),gps(iy),gps(iz),orientation_list(iorientation,:),[0 0 0],eX,eY,eZ);
                test(iorientation,ix,iy,iz,:)=get_f_dervs_3d(gps(ix),gps(iy),gps(iz),orientation_list(iorientation,:),[0 0 0]);
                
                %----
                for i=2:1:6
                    diff=abs(test(iorientation,ix,iy,iz,i)-weights_plus(iorientation,1,i+1,ix,iy,iz));
                    if diff>10^-10
                        iorientation
                        [ix iy iz]
                        i
                        test(iorientation,ix,iy,iz,i)
                        weights_plus(iorientation,1,i+1,ix,iy,iz)
                        test(iorientation,ix,iy,iz,:);
                        error('found')
                    end
                end
                %---
                
            end
        end
    end
end

%--------------------------------------------- type 000, order 000
for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                sol(iorientation,1,1,ix,iy,iz)=weights_plus(iorientation,1,1,ix,iy,iz);
            end
        end
    end
end
%--------------------------------------------- type 000, order~000
for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                
%                 %---
%                 for iorder=2:1:7
%                     val=0;
%                     irow=iorder-1;
%                     for j=1:1:9
% %                         val=val+invJ(ix,iy,iz,irow,j)*weights_plus(iorientation,1,j+1,ix,iy,iz);
%                         val=val+invJ(ix,iy,iz,irow,j)*test(iorientation,ix,iy,iz,irow);
%                     end
%                     %---------------------------------------------
%                     if sol(iorientation,1,iorder,ix,iy,iz)~=0
%                         error('Over computed')
%                     end
%                     %---------------------------------------------
%                     sol(iorientation,1,iorder,ix,iy,iz)=val;
% %                     get_fg_dervs_3d(alpha,beta,gamma,orientation_list(iorientation,:),types,eX,eY,eZ)
%                 end
%                 %---
                
                gf_dervs=get_fg_dervs_3d(gps(ix),gps(iy),gps(iz),orientation_list(iorientation,:),[0 0 0],eX,eY,eZ);
                
                for i=1:1:6
                    sol(iorientation,1,i+1,ix,iy,iz)=gf_dervs(i);
                end
                
            end
        end
    end
end

%---------------------------------------------------------------------
%-----UPDATED TODAY--------------------------------------------
%--------------------------------------------- type 100, order 000
% return
% map_to_global_1_compt_3d(alpha,beta,gamma,eCompt,orders)

derivative_terms=derivative_terms_3d(gps,eX,eY,eZ);

for iorientation=1:1:8
    for itype=2:1:4
        icompt=itype-1;
        for iz=1:1:3
            for iy=1:1:3
                for ix=1:1:3
                    sol(iorientation,itype,1,ix,iy,iz)=...
                        weights_plus(iorientation,itype,1,ix,iy,iz)...
                        *derivative_terms(icompt,ix,iy,iz);
                end
            end
        end
    end
end

%--------------------------------------------- type 110, order 000

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                sol(iorientation,5,1,ix,iy,iz)=...
                    double_Gaussian_quadrature(@eval_integrant_type5_order0,...
                    gps(ix),gps(iy),gps(iz),gps,w,orientation_list(iorientation,:),eX,eY);
                sol(iorientation,6,1,ix,iy,iz)=...
                    double_Gaussian_quadrature(@eval_integrant_type6_order0,...
                    gps(ix),gps(iz),gps(iy),gps,w,orientation_list(iorientation,:),eX,eZ);
                 sol(iorientation,7,1,ix,iy,iz)=...
                    double_Gaussian_quadrature(@eval_integrant_type7_order0,...
                    gps(iy),gps(iz),gps(ix),gps,w,orientation_list(iorientation,:),eY,eZ);
            end
        end
    end
end

%--------------------------------------------- type 100, order ~000

gf_dervs=get_fg_dervs_3d(gps(1),gps(1),gps(1),[0 0 0],[1 0 0],eX,eY,eZ);
for i=1:1:6
    sol(1,2,i+1,1,1,1)=gf_dervs(i);
end

%-------------------- The following is corrct
sol(1,2,2,1,1,1)=local_basis_3d(gps(1),gps(1),gps(1),[0 0 0],[1 0 0],[1 0 0]);

% 
% get_invJ(gps(ix),gps(iy),gps(iz),eX,eY,eZ);

% get_fg_dervs_3d(alpha,beta,gamma,orientations,types,eX,eY,eZ);
% dlmwrite('test_file',eZ);

% sol(1,5,2,1,1,1)=double_Gaussian_quadrature(@test198321,gps(1),gps(1),gps(1),gps,w,[0 0 0],eX,eY);


%=====================================

% alpha
% beta
% gamma

% gps1=gps*gps(1);
% gps2=gps1;
% gamma=gps(1);
% 
% summation=0;
% 
% for i=1:1:3
%     for j=1:1:3
%         
%         alpha=gps1(i);
%         beta=gps2(j);
%         
%         term1=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
%         term2=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
%         term3=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 1 0]);
%         term4=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[1 0 0]);
%         
%         det=term1*term2-term3*term4;
%         
%         take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 1 0],eX,eY,eZ);
%         
%         summation=summation+w(i)*w(j)*take(1)*det;
%         
%     end
% end
% summation=summation*gps(1)^2;

%--------------- GOOD ONE---------------------
gps1=gps*gps(1);

summation=0;
beta=gps(1);
gamma=gps(1);
for i=1:1:3
    alpha=gps1(i);
    take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 0 0],eX,eY,eZ);
    summation=summation+w(i)*take(3);
end

sol(1,2,4,1,1,1)=summation*map_to_global_1_compt_3d(gps(1),gps(1),gps(1),eX,[1 0 0])*gps(1);
% sol(1,2,4,1,1,2)=map_to_global_1_compt_3d(gps(1),gps(1),gps(1),eX,[1 0 0]);
%-----------------------------------------------------------

%--------------- GOOD ONE---------------------
gps1=gps*gps(1);
summation=0;
beta=gps(1);
gamma=gps(1);
for i=1:1:3
    alpha=gps1(i);
    take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 0 0],eX,eY,eZ);
    summation=summation+w(i)*take(6);
end
sol(1,2,7,1,1,1)=summation*map_to_global_1_compt_3d(gps(1),gps(1),gps(1),eX,[1 0 0])*gps(1);
%-----------------------------------------------------------

%--------------- GOOD ONE---------------------
gps1=gps*gps(1);
summation=0;
alpha=gps(1);
gamma=gps(1);
for i=1:1:3
    beta=gps1(i);
%     take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 0 0],eX,eY,eZ);
    tale=local_basis_3d(alpha,beta,gamma,[0 0 0],[1 1 0],[1 1 0]);

    summation=summation+w(i)*tale;
end
sol(1,5,2,1,1,1)=summation*map_to_global_1_compt_3d(gps(1),gps(1),gps(1),eY,[0 1 0])*gps(1);

if sol(1,5,2,1,1,1)==0
    error('adfadsf')
end
%-----------------------------------------------------------

%--------------- GOOD ONE---------------------

gps1=gps*gps(1);
gps2=gps*gps(1);
gamma=gps(1);
summation=0;
itest=0;
test101=zeros(1,9);

for i=1:1:3
    for j=1:1:3
        
        itest=itest+1;
        
        alpha=gps1(i);
        beta=gps2(j);
        
        term1=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[1 0 0]);
        
        det=term1*term2-term3*term4;
        
        if det==0
            error('sfadfs')
        end
        
        take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 1 0],eX,eY,eZ);
        
        summation=summation+w(i)*w(j)*take(3)*det;
        
        test101(itest)=summation;
        
    end
end
summation=summation*gps(1)^2;

sol(1,5,4,1,1,1)=summation;

%-----------------------------------------------------------


gps1=gps*gps(1);
summation=0;
alpha=gps(1);
gamma=gps(1);
for i=1:1:3
    beta=gps1(i);
%     take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 0 0],eX,eY,eZ);
%     tale=local_basis_3d(alpha,beta,gamma,[0 0 0],[1 1 0],[1 1 0]);
    
    take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 1 0],eX,eY,eZ);

    summation=summation+w(i)*take(1);
end
sol(1,5,5,1,1,1)=summation*map_to_global_1_compt_3d(gps(1),gps(1),gps(1),eY,[0 1 0])*gps(1);



test101

warning('adfasfd')

% local_basis_3d(alpha,beta,gamma,orientations,types,orders)

% teafda=sol(1,2,4,1,1,1,2)
% sol(1,5,2,1,1,1)=summation;

take=get_fg_dervs_3d(gps(1),gps(1),gps(1),[0 0 0],[1 0 0],eX,eY,eZ);

sol(1,2,5,1,1,1)=take(1);

end

%=====================================
%=====================================
%=====================================
%=====================================
%=====================================
%=====================================
%=====================================
%=====================================
%=====================================

function sol=test198321(alpha,beta,gamma,orientations,eX,eY)

eZ=dlmread('test_file');

take=get_fg_dervs_3d(alpha,beta,gamma,orientations,[1 1 0],eX,eY,eZ);

term1=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
term2=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
term3=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 1 0]);
term4=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[1 0 0]);

det=term1*term2-term3*term4;
sol=take(1)*det;

end

function test234245



end

%=====================================
%=====================================

function sol=double_Gaussian_quadrature(fun,val1,val2,val3,gps,w,orientations,compt1,compt2)
% It evaluates from 0 to val1 on one component and 0 to val2 on another component
gps1=gps*val1;
gps2=gps*val2;
sol=0;
for i=1:1:3
    for j=1:1:3
        sol=sol+w(i)*w(j)*fun(gps1(i),gps2(j),val3,orientations,compt1,compt2);
    end
end
sol=sol*val1*val2;
end


% function sol=test_integrant(alpha,beta,gamma,eX,eY,eZ)
% term1=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
% term2=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
% term3=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 1 0]);
% term4=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[1 0 0]);
% B_derv=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[1 1 0],eX,eY,eZ);
% sol=B_derv*(term1*term2-term3*term4);
% end

function sol=eval_integrant_type7_order0(beta,gamma,alpha,orientations,eY,eZ)
term1=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
term2=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,[0 0 1]);
term3=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 0 1]);
term4=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,[0 1 0]);
B=local_basis_3d(alpha,beta,gamma,orientations,[0 1 1],[0 1 1]);
sol=B*(term1*term2-term3*term4);
end

function sol=eval_integrant_type6_order0(alpha,gamma,beta,orientations,eX,eZ)
term1=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
term2=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,[0 0 1]);
term3=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 0 1]);
term4=map_to_global_1_compt_3d(alpha,beta,gamma,eZ,[1 0 0]);
B=local_basis_3d(alpha,beta,gamma,orientations,[1 0 1],[1 0 1]);
sol=B*(term1*term2-term3*term4);
end

function sol=eval_integrant_type5_order0(alpha,beta,gamma,orientations,eX,eY)
term1=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
term2=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
term3=map_to_global_1_compt_3d(alpha,beta,gamma,eX,[0 1 0]);
term4=map_to_global_1_compt_3d(alpha,beta,gamma,eY,[1 0 0]);
B=local_basis_3d(alpha,beta,gamma,orientations,[1 1 0],[1 1 0]);
sol=B*(term1*term2-term3*term4);
end


%=====================================
%=====================================
%=====================================

function sol=get_fg_dervs_3d(alpha,beta,gamma,orientations,types,eX,eY,eZ)
invJ=get_invJ(alpha,beta,gamma,eX,eY,eZ);
f_der=get_f_dervs_3d(alpha,beta,gamma,orientations,types);
sol=zeros(1,6);
for i=1:1:6
    val=0;
    for j=1:1:9
        val=val+invJ(i,j)*f_der(j);
    end
    sol(i)=val;
end

end

function f_der=get_f_dervs_3d(alpha,beta,gamma,orientations,types)
f_der=zeros(1,9);
f_der(1)=local_basis_3d(alpha,beta,gamma,orientations,types,[1 0 0]+types);
f_der(2)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 1 0]+types);
f_der(3)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 0 1]+types);
f_der(4)=local_basis_3d(alpha,beta,gamma,orientations,types,[2 0 0]+types);
f_der(5)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 2 0]+types);
f_der(6)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 0 2]+types);
f_der(7)=local_basis_3d(alpha,beta,gamma,orientations,types,[1 1 0]+types);
f_der(8)=local_basis_3d(alpha,beta,gamma,orientations,types,[1 0 1]+types);
f_der(9)=local_basis_3d(alpha,beta,gamma,orientations,types,[0 1 1]+types);
end

function sol=derivative_terms_3d(gps,eX,eY,eZ)
sol=zeros(3,3,3,3);
for iz=1:1:3
    for iy=1:1:3
        for ix=1:1:3
            sol(1,ix,iy,iz)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 0 0]);
            sol(2,ix,iy,iz)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
            sol(3,ix,iy,iz)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[0 0 1]);
        end
    end
end
end

function sol=local_basis_3d(alpha,beta,gamma,orientations,types,orders)
sol=basis(alpha,orientations(1),types(1),orders(1))...
    *basis(beta,orientations(2),types(2),orders(2))...
    *basis(gamma,orientations(3),types(3),orders(3));
end

function invJ=e_invJ(gps,eX,eY,eZ)
invJ=zeros(3,3,3,9,9);
for iz=1:1:3
    for iy=1:1:3
        for ix=1:1:3
            invJ(ix,iy,iz,:,:)=get_invJ(gps(ix),gps(iy),gps(iz),eX,eY,eZ);
        end
    end
end
end