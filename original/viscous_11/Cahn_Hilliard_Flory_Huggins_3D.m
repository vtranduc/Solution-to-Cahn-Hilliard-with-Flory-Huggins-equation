function Cahn_Hilliard_Flory_Huggins_3D

clear
clc

tic
counter_init=cputime;

% ===================== User defined variables ============================

nex=2;
ney=1;
nez=1;

dto=2.0e-6;
dt=2.0e-6;

obs_t=2.0e-1; %Observation time

diff=1800;

T=0.6;

ci = 0.7; %Overall initial concentrations
ci_fluc = 0.01; %Initial fluctuation

nr_tol=1.0e-6;
nr_max_iteration=20;

n1=1;
n2=10;

time_out=3600*6;

entropy=1;
T_theta=1;

fps=10;

dt_down=0.5;
dt_up=1.2;
dt_ideal=dt;
bypass_tol=5;

max_frame=100*100;

% xlen=nex;
% ylen=ney;
% zlen=nez;

xlen=1;
ylen=1;
zlen=1;

% === Specify constants which will be used throughout =====================

w=[0.27778 0.4444 0.27778];
gp=[0.1127 0.5 0.8873];

nx=nex+1;ny=ney+1;nz=nez+1;
ne=nex*ney*nez;
n=nx*ny*nz;
neight=8*n;

ll=2;
li=8;
lu=2+8*ney;
bl=3;
bi=8*ny;
bu=3+8*nex*ny;
tl=3+8*ney;
ti=8*ny;
tu=3+8*ney+8*nex*ny;
rl=2+8*nex*ny;
ri=8;
ru=2+8*ney+8*nex*ny;

zi=nx*ny*8;zu=zi*(nz-1);bny=nx*ny*8

% === initialize time =====================================================

time=0.0;

% === Evaluate value of chi ===============================================

chi=0.5-entropy*(1-T_theta/T);
% ======================= Start the counter ===============================


% === Determine the coordinate at each node ===============================

x=zeros(1,n);y=zeros(1,n);z=zeros(1,n);

x_space=xlen/nex;
for i=1:nx
    x_coord=x_space*(i-1);
    for j=1+ny*(i-1):nx*ny:1+ny*(i-1)+nx*ny*(nz-1)
        x(j:j+ny-1)=x_coord;
    end
end

y_space=ylen/ney;
for i=1:ny
    y_coord=y_space*(i-1);
    for j=1:nz
        for k=i+nx*ny*(j-1):ny:i+nx*ny*(j-1)+ny*(nx-1)
            y(k)=y_coord;
        end
    end
end

z_space=zlen/nez;
for i=1:nz
    z_coord=z_space*(i-1);
    z(nx*ny*(i-1)+1:nx*ny*(i-1)+nx*ny)=z_coord;
end

% === Generate matrices that help to characterize each element uniquely ===

% nop=zeros(ne,16);
% for i=1:ne
%     nop(i,1)=4*floor((i-1)/ney)+4*(i-1)+1;
%     nop(i,9)=nop(i,1)+4*ny;
%     for j=1:7
%         nop(i,1+j)=nop(i,1)+j;
%         nop(i,9+j)=nop(i,9)+j;
%     end
% end
% nopm=zeros(i,4);
% for i=1:ne
%     nopm(i,1)=floor((i-1)/ney)+i;
%     nopm(i,2)=nopm(i,1)+1;
%     nopm(i,3)=nopm(i,1)+ny;
%     nopm(i,4)=nopm(i,3)+1;
% end

% === Randomize and apply initial conditions ==============================

co = zeros(1,neight);

for i=1:neight
    co(i)=ci+(ci_fluc*(2*rand(1,1)-1));
end

co=icc3D(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny);

% itest=16
% co(itest*8-7:itest*8)'
%%{
% === Set up Residual vector and Jacobian matrix ==========================

% === Take the very first time step =======================================

c=co;

% === Set up for plotting =================================================

% === Plot the concentration in second step, which is qual to first one ===

% === Prepare to enter main loop ==========================================

cp=zeros(1,neight);
nodes=zeros(1,8);
gbasis=zeros(1,64);

nxy=nx*ny;

% === Enter timer loop ====================================================

for rahmat=[1]
    
% === Update the concentration data =======================================
    
    for i=1:1:neight
        cp(i)=c(i)+dt*((c(i)-co(i))/dto);
    end
    coo=co;
    co=c;
    c=cp;
    
% === Apply BC to newly predicted concentration ===========================

    c=bc3D(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny);
    
% === Enter Newton-Raphson iterations =====================================
    
    for saeedi=1:1:4
        
        sf=zeros(1,neight);
        sj=zeros(neight,neight);
        
% === Start looping through each element ==================================
        
        for e=1:1:ne

% === Characterize the current element ====================================
            
            z1=mod(e,ney);
            if z1==0
                z1=ney;
            end
            z2=(e-z1)/ney;
            z3=mod(z2,nex);
            z4=(z2-z3)/nex;
            
            nodes(1)=nx*ny*z4+ny*z3+z1;
            nodes(2)=nodes(1)+1;
            nodes(3)=nodes(1)+ny;
            nodes(4)=nodes(3)+1;
            nodes(5)=nodes(1)+nxy;
            nodes(6)=nodes(5)+1;
            nodes(7)=nodes(5)+ny;
            nodes(8)=nodes(7)+1;
            
            for i=1:1:8
                z1=(i-1)*8+1;
                gbasis(z1:z1+7)=nodes(i)*8-7:1:nodes(i)*8;
            end
            
            dx=x(nodes(3))-x(nodes(1));
            dy=y(nodes(2))-y(nodes(1));
            dz=z(nodes(5))-z(nodes(1));
            
% === Shift the weights of two basis functions ============================
            
            for w1=1:1:3
                for w2=1:1:3
                    for w3=1:1:3

% === Characterize the current element's basis accordingly ================

                        weights=...
                            tfunct3D([gp(w1),gp(w2),gp(w3)],[dx dy dz]);
                    
                        %YOU MUST RESOLVE WHAT TO DO WITH THE
                        %TRANSFORMATION!!!!
                        
% === Determine previous and current absolute values, and associated slopes
% at current particular pint ==============================================

                        con=0.0;
                        cono=0.0;
                        conx=0.0;
                        cony=0.0;
                        conz=0.0;
                        conxx=0.0;
                        conyy=0.0;
                        conzz=0.0;
                        
                        for i=1:1:64
                            con=con+c(gbasis(i))*weights(i,1);
                            cono=cono+co(gbasis(i))*weights(i,1);
                            conx=conx+c(gbasis(i))*weights(i,2);
                            cony=cony+c(gbasis(i))*weights(i,3);
                            conz=conz+c(gbasis(i))*weights(i,4);
                            conxx=conxx+c(gbasis(i))*weights(i,5);
                            conyy=conyy+c(gbasis(i))*weights(i,6);
                            conzz=conzz+c(gbasis(i))*weights(i,7);
                        end
                        
% === Change in concentration over time step is estimated =================

                        cont=(con-cono)/dt;
                        
% === Fill in residual vector =============================================

                        sq_term=conx^2+cony^2+conz^2;
                        cLaplace=conxx+conyy+conzz;
                        first_dev_sum=conx+cony+conz;
                        chi_term=(1/(con*n1)+1/((1-con)*n2)-2*chi);
                        
                        for i=1:1:64
                            
                            iphiLaplace=0.0;
                            for j=5:1:7
                                iphiLaplace=iphiLaplace+weights(i,j);
                            end
                            
                            iphi=weights(i,1);
                            
                            sf(gbasis(i))=sf(gbasis(i))-w(w1)*w(w2)*w(w3)*dx*dy*dz*...
                                (iphi*cont-diff*T*iphi*((-1/(con^2*n1)+1/((1-con)^2*n2))*...
                                sq_term+chi_term*cLaplace)+...
                                cLaplace*iphiLaplace);
                            
                            for j=1:1:64
                                
                                jphi=weights(j,1);
                                
                                jphiLaplace=0.0;
                                for k=5:1:7
                                    jphiLaplace=jphiLaplace+weights(j,k);
                                end
                                sj(gbasis(i),gbasis(j))=sj(gbasis(i),gbasis(j))+...
                                    w(w1)*w(w2)*w(w3)*dx*dy*dz*...
                                    (iphi*jphi/dt-diff*T*iphi*((2*jphi/(con^3*n1)+...
                                    2*jphi/((1-con)^3*n2))*sq_term+...
                                    (-1/(con^2*n1)+1/((1-con)^2*n2))*2*jphi*first_dev_sum+...
                                    (-jphi/(con^2*n1)+jphi/((1-con)^2*n2))*cLaplace+...
                                    chi_term*jphiLaplace)+...
                                    jphiLaplace*iphiLaplace);
                            end
                        end
                    end
                end
            end
        end
        
        % === Apply boundary conditions to sf and sj ==============================
        fprintf('Apply bc')
        [sf,sj]=bsfsfsj3D(sf,sj,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny,neight);
        
        % === Create compact matrix ===============================================
        
%         fprintf('Compressing')
%         toc
%         %[rowi,coli,val]=compress_sj(sj,nfour,nfour);
%         [rowi,coli,val]=find(sj);
%         toc
%         %[rowi',coli',val']
%         fprintf('Sparsing')
%         sj=sparse(rowi,coli,val);
%         toc
        
        % === Carry matrix Division ===============================================
        determinant=det(sj)
%         if determinant==0
%             dlmwrite('seeMatrix.xlsx',sj)
%             itry=0;
%             for iii=64:1:neight
%                 if prod(get_diagonal(iii,sj))~=0
%                     iii
%                     z1=get_diagonal(iii,sj);
%                     z1'
%                     error('Non zero found!')
%                 end
%             end
%             itry
%             error('SINGULAR!!!')
%         end
        fprintf('Carry out division')
        c_= sj \ sf';
        
        % === Update the solution =================================================

        c = c + c_';
        
% === Evaluate the error ==================================================
        fprintf('error eval')
        err = sqrt(sum(c_.^2.0))
        
    end
end

sf';
sf=icc3D(sf,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny);
sf';

% itest=36;
% jtest=31;
% 
% sj(itest*8-6,jtest*8-7)

%%}

end