function test58

clear
clc


% === Start the timer =====================================================

counter_init=tic;

% === User defined variables ==============================================

% --- Define mesh ---------------------------------------------------------

nex=5;
ney=4;
nez=3;

shape=3;

% So far, we have the following available
% 1 - Cube
% 2- Sphere

% The data must have been made ready with the following file

% main_isoparametric_mapping_3d

% --- Polymer properties --------------------------------------------------

n1=1;
n2=10;
 
diff=5000;

%======================
ci=0.45;
T=[0.4];
%======================

include_fluc=1;
ci_fluc = 0.0001;

entropy=1;
T_theta=1;

ci_critical=0; % This bypasses ci specified above
T_critical=0; % This bypasses T specified above

% --- Simulation properties -----------------------------------------------

parallel_computing=1; % Non parallel computing uses traditional method

show_figure=0;
export_data=0; % Yes=1, No=0


dt=1.0*10^-7;

obs_t=2.0e-1;

nr_tol=1.0e-6;
nr_max_iteration=20;
%======================
time_out=3600*0.3;
max_frame=1000;
%======================
fps=10;
dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-150;
dt_ideal=5.0e-10;
bypass_tol=5;

tol=1.0e-6;
sparse_traditional_sj=0;

predict_with_acceleration=1;

commit_sfsj_terms_to_RAM=1;

% --- Phase diagram properties --------------------------------------------

nT=100;
nInitialGuessPts=100;
guess_accuracy=10^-3;

proceed_request=1;

ifig_phase=1;

% --- Analysis specifications ---------------------------------------------

alpha=200;
beta=8;

% --- Illustration specifications -----------------------------------------

frame_pos=[0 0];
frame_size=[720 720];
ifig_simulation=2;

face_colors={'yellow','blue','red'};
% The order is {Mean Lower Higher}
view_angle=3;
transparency=0.9;
edge_color='none';
ifig_conclusion=3;

domain_EdgeColor='green';
domain_FaceColor='white';
domain_face_transparency=0.1;
domain_edge_transparency=0.8;

% === Validate inputs =====================================================

% VALIDATE THEM HERE

if shape==3
    if parallel_computing~=1 || commit_sfsj_terms_to_RAM~=1
        error('Only one method is available for shape 3')
    end
end
% return

% === Set up variables ====================================================

[nx,ny,nz,n,neight]=setUp_3d(nex,ney,nez);
if ci_critical==1 && T_critical==1
    [ci,T]=identify_critical_point(n1,n2,entropy,T_theta);
elseif ci_critical==1
    [ci,~]=identify_critical_point(n1,n2,entropy,T_theta);
elseif T_critical==1
    [~,T]=identify_critical_point(n1,n2,entropy,T_theta);
end

c_obs=[mean(ci) mean(ci)-ci_fluc mean(ci)+ci_fluc];
ne=nex*ney*nez;
nexney=nex*ney;
dto=inf;
if predict_with_acceleration==1
    dtoo=inf;
end
bypass=0;



if shape==1 || shape==2
    nnz_=nnz_sj_3d(nx,ny,nz);
    [gbfs1,gbfs2,iZeros,iOnes]=sj_mapper_3d(nx,ny,nz);
    fprintf('Importing Data...\n')
    [X,Y,Z,weights]=...
        retrieve_isoparametrically_mapped_weights_3d(...
        shape,nex,ney,nez,parallel_computing);
    fprintf('Computing determinants...\n')
    determinants=isoparametric_determinant_3d(...
        ne,nex,nx,ny,nexney,X,Y,Z,parallel_computing);
elseif shape==3
    nxny=nx*ny;
    nx_=nx-2;
    ny_=ny-2;
    [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
    neTotal=ne+n_eSurf(6);
    nTotal=n+n_nSurf(6);
    nnz_=nnz_sj_sph(nx,ny,nz,n_nSurf);
    [gbfs1,gbfs2]=sj_mapper_sph(n,nx,ny,nz,nnz_,nxny,n_nSurf,nx_,nTotal);
    [n_nSurf,n_eSurf]=get_n_nSurf_n_eSurf(nex,ney,nez,ny,nz,nx_,ny_);
    fprintf('Importing Data...\n')
    [X,Y,Z,weights]=...
        retrieve_isoparametrically_mapped_weights_3d(...
        shape,nex,ney,nez,parallel_computing);
    fprintf('Computing determinants...\n')
    
    %-----Test---------------
    
%     size(weights)
%     
%     return

    take=dlmread('test_weights_surface.dat');
    weights_surface=zeros(n_eSurf(6),8,8,7,3,3);
    
    disp('reformatting')
    toc(counter_init)
    
    for e=1:1:n_eSurf(6)
        index=0;
        eInfo=take(e,:);
        for iorientation=1:1:8
            for itype=1:1:8
                for iorder=1:1:7
                    for idim1=1:1:3
                        for idim2=1:1:3
                            index=index+1;
                            weights_surface(e,iorientation,itype,iorder,idim1,idim2)=...
                                eInfo(index);
                        end
                    end
                end
            end
        end
    end
    
    
    wTerms_surface=generate_wTerms_curved_surface_sph(n_eSurf,weights_surface);
    
    determinants_surface=isoparametric_determinant_surface_sph(...
        ne,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,n_eSurf,n_nSurf,...
        nx_,parallel_computing);
    
    take=determinants_surface(:,1,1);
    nnz(take)
    size(take)
    n_eSurf(6)
    
    
    
    

    disp('----new test--------')
%     return
    
    normals=get_surface_normal_sph(ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_,X,Y,Z,parallel_computing);
    
%     normals(1,:)
%     
%     hold on
%     for i=1:1:n_eSurf(6)
%         plot3([normals(i,1),0],[normals(i,2),0],[normals(i,3),0])
%         
%     end
%     hold off
%     grid on
%     xlabel('x');ylabel('y');zlabel('z')
%     
%     return
    %--------------------------
    
    
    determinants=isoparametric_determinant_sph(...
        ne,neTotal,nex,ney,n,nx,ny,nz,nxny,nexney,X,Y,Z,...
        n_eSurf,n_nSurf,nx_,parallel_computing);
    nTotal=n+n_nSurf(6);
    neightTotal=neight+8*n_nSurf(6);
end


toc(counter_init)
fig_domain=[min(X) max(X) min(Y) max(Y) min(Z) max(Z)];
% THIS FUNCTION CAN BE PUT IN DIFFERENT MAIN FILE
toc(counter_init)
if parallel_computing==1
    if shape==1 || shape==2
        wTerms=generate_wTerms_curved_3d(ne,weights);
    elseif shape==3
        wTerms=generate_wTerms_curved_3d(neTotal,weights);
    end
end


if show_figure==1
    [x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot]=...
        setUp_illustration_curved_structure_3d(nx,ny,nz,X,Y,Z);
end

% --- Break down temperature gradient -------------------------------------

if shape==1 || shape==2
    [diffT,chi]=create_temperature_gradient_3d(...
        parallel_computing,shape,T,diff,entropy,T_theta,nx,ny,nexney,ne,nex,X,Y,Z);
elseif shape==3
    [diffT,chi]=create_temperature_gradient_3d(...
        parallel_computing,shape,T,diff,entropy,T_theta,nx,ny,nexney,neTotal,nex,X,Y,Z);
end

% --- Set up folder for outputting results --------------------------------

if export_data==1 || show_figure==1
    output_path=folder_setUp_3d(export_data);
    
    if export_data==1
        sf=[nex;ney;nez;dt;obs_t;diff;include_fluc;ci_fluc;nr_tol;...
            nr_max_iteration;time_out;max_frame;entropy;T_theta;...
            dt_down;dt_up;dt_min;dt_ideal;bypass_tol;n1;n2];
        dlmwrite([output_path 'specification/spec_general'],sf);
        dlmwrite([output_path 'specification/status'],0);
        dlmwrite([output_path 'specification/spec_T'],T);
        dlmwrite([output_path 'specification/spec_ci'],ci);
    end

% === Phase diagram =======================================================

    fprintf('Generating phase diagram...\n')
    [xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
        n1,n2,entropy,T_theta,nT,tol,parallel_computing,nInitialGuessPts,...
        guess_accuracy);
	if show_figure==1
        figure(ifig_phase)
        illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
        if proceed_request==1
            proceed=proceed_dialog();
            if proceed~=1
                return
            end
        end
	end
end

% === Generate initial condition ==========================================

%THIS FUNCTION IS NOT DESGINED FOR RADIAL GRADIENT!!!   
if shape==1 || shape==2
    co=generate_co_3d(ci,include_fluc,ci_fluc,n,neight,nx,ny,nz);
elseif shape==3
    co=generate_co_3d(ci,include_fluc,ci_fluc,nTotal,neightTotal,nx,ny,nz);
end

% === Initialize the first time step ======================================

time=0.0;
c=co;

if predict_with_acceleration==1
    coo=co;
end

iframe=1;

% --- Plot the first time step --- %

if show_figure==1
    fig=figure(ifig_simulation);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    clf(fig)
    visual(iframe)=illustrate_and_analyze_curved_structure_3d(...
        c,c_obs,fig_domain,face_colors,edge_color,...
        nx,ny,nz,view_angle,transparency,time,...
        x_minus_x,x_minus_y,x_minus_z,x_plus_x,x_plus_y,x_plus_z,...
        y_minus_x,y_minus_y,y_minus_z,y_plus_x,y_plus_y,y_plus_z,...
        z_minus_x,z_minus_y,z_minus_z,z_plus_x,z_plus_y,z_plus_z,...
        x_plot,y_plot,z_plot,domain_EdgeColor,domain_FaceColor,...
        domain_face_transparency,domain_edge_transparency);
end

% --- Output the first time step --- %

if export_data==1
    dlmwrite([output_path 'time/iteration_1'],0);
    dlmwrite([output_path 'concentration/iteration_1'],c);
end

% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================
% ==== Enter main loop ====================================================


fprintf('Enter main cycle\n')
cont=zeros(3,3,3);
w=[5/18 4/9 5/18];

for fdasfasdfa=1:1:5
    
    sf=zeros(neightTotal,1);
    sj=zeros(neightTotal,neightTotal);
    
    
    conc_=get_conc_curved_sph(neTotal,c,weights,ne,nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);
    cono=get_conc_type_I_curved_sph(neTotal,co,weights,ne,...
        nexney,nex,ney,nx,ny,nz,nxny,n,n_eSurf,n_nSurf,nx_);

    
%     [c' co']
%     
%     size(conc_)
%     
%     take=conc_(:,1,:,:,:);
%     
%     
%     mean(mean(mean(mean(mean(mean(mean(mean(mean(mean(take))))))))))
%     mean(mean(mean(mean(mean(mean(mean(mean(mean(mean(cono))))))))))
    
%     return

    %===============================================================
%     for e=1:1:neTotal
    for e=116
        
        
        [e,ne+n_eSurf(6)]
        gbfs=get_gbfs_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
        
        
        for i=1:1:64
            [node,type]=analyze_gbs_3d(gbfs(i));
            fprintf('%i %i %i\n',i,node,type)
        end
        
        size(weights)
        
        xp=1;
        yp=2;
        zp=2;
        
        orderp=7;
        
        take=0;
        iorientation=1;
        itype=0;
        for i=1:1:64
            itype=itype+1;
            if itype==9
                itype=1;
                iorientation=iorientation+1;
            end
            [iorientation itype]
            take=take+c(gbfs(i))*weights(e,iorientation,itype,orderp,xp,yp,zp);
        end
        
        take
        
        conc_(e,orderp,xp,yp,zp)
        
        size(conc_)
        
        size(cono)
        
        cono(e,xp,yp,zp)
        conc_(e,1,xp,yp,zp)
        
        
        nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
        hold on
        for i=1:1:8
%             [X(nodes(i)),Y(nodes(i)),Z(nodes(i))]
            plot3(X(nodes(i)),Y(nodes(i)),Z(nodes(i)),'ro')
        end
        grid on
        hold off
        
        
        
        weights(e,1,1,1,1,1,1)
        weights(e,1,1,1,2,2,2)
        weights(e,1,1,1,3,3,3)
        
        return
        
        for iz=1:1:3
            for iy=1:1:3
                for ix=1:1:3
                    cont(ix,iy,iz)=(conc_(e,1,ix,iy,iz)-cono(e,ix,iy,iz))/dt;
                end
            end
        end
        for iz=1:1:3
            for iy=1:1:3
                for ix=1:1:3
                    isf=0;
                    for iorientation=1:1:8
                        for itype=1:1:8
                            isf=isf+1;
                            sf(gbfs(isf))=sf(gbfs(isf))+determinants(e,ix,iy,iz)*w(ix)*w(iy)*w(iz)*(...
                                weights(e,iorientation,itype,1,ix,iy,iz)*cont(ix,iy,iz)-diffT*weights(e,iorientation,itype,1,ix,iy,iz)*...
                                ((-1/(conc_(e,1,ix,iy,iz)^2*n1)+1/((1-conc_(e,1,ix,iy,iz))^2*n2))*...
                                (conc_(2,ix,iy,iz)^2+conc_(3,ix,iy,iz)^2+conc_(4,ix,iy,iz)^2)+...
                                (1/(conc_(e,1,ix,iy,iz)*n1)+1/((1-conc_(e,1,ix,iy,iz))*n2)-2*chi)*...
                                (conc_(5,ix,iy,iz)+conc_(6,ix,iy,iz)+conc_(7,ix,iy,iz)))+...
                                (conc_(5,ix,iy,iz)+conc_(6,ix,iy,iz)+conc_(7,ix,iy,iz))*...
                                (weights(e,iorientation,itype,5,ix,iy,iz)+weights(e,iorientation,itype,6,ix,iy,iz)+weights(e,iorientation,itype,7,ix,iy,iz))...
                                );
                            jsj=0;
                            %--------------------
                            for jorientation=1:1:8
                                for jtype=1:1:8
                                    jsj=jsj+1;

                                    sj(gbfs(isf),gbfs(jsj))=sj(gbfs(isf),gbfs(jsj))+determinants(e,ix,iy,iz)*w(ix)*w(iy)*w(iz)*(...
                                        weights(e,iorientation,itype,1,ix,iy,iz)*weights(e,jorientation,jtype,1,ix,iy,iz)/dt...
                                        -diffT*weights(e,iorientation,itype,1,ix,iy,iz)...
                                        *(2*weights(e,jorientation,jtype,1,ix,iy,iz)...
                                        *((conc_(e,1,ix,iy,iz)^3*n1)^-1+((1-conc_(e,1,ix,iy,iz))^3*n2)^-1)...
                                        *(conc_(2,ix,iy,iz)^2+conc_(3,ix,iy,iz)^2+conc_(4,ix,iy,iz)^2)...
                                        +weights(e,jorientation,jtype,1,ix,iy,iz)...
                                        *(-(conc_(e,1,ix,iy,iz)^2*n1)^-1+((1-conc_(e,1,ix,iy,iz))^2*n2)^-1)...
                                        *(2*(conc_(2,ix,iy,iz)+conc_(3,ix,iy,iz)+conc_(4,ix,iy,iz))...
                                        +(conc_(5,ix,iy,iz)+conc_(6,ix,iy,iz)+conc_(7,ix,iy,iz)))...
                                        +((conc_(e,1,ix,iy,iz)*n1)^-1+((1-conc_(e,1,ix,iy,iz))*n2)^-1-2*chi)...
                                        *(weights(e,jorientation,jtype,5,ix,iy,iz)+weights(e,jorientation,jtype,6,ix,iy,iz)+weights(e,jorientation,jtype,7,ix,iy,iz)))...
                                        +(weights(e,iorientation,itype,5,ix,iy,iz)+weights(e,iorientation,itype,6,ix,iy,iz)+weights(e,iorientation,itype,7,ix,iy,iz))...
                                        *(weights(e,jorientation,jtype,5,ix,iy,iz)+weights(e,jorientation,jtype,6,ix,iy,iz)+weights(e,jorientation,jtype,7,ix,iy,iz))...
                                        ...
                                        );
                                end
                            end
                            %------------
                        end
                    end
                end
            end
        end
    end
    %==========================================================================================
    
    for ijk=1:1:n_nSurf(6)
%         how=how+1;
        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+2,neightTotal);
        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+5,neightTotal);
        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+6,neightTotal);
        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+8,neightTotal);

        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+3,neightTotal);
        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+4,neightTotal);
        [sf,sj]=raw_boundary_applier(sf,sj,n*8+(ijk-1)*8+7,neightTotal);
    end
    
    sf'
    
    
    c_=sj\-sf;
    c=c+c_';
    err=sqrt(sum(c_.^2.0))

end









end