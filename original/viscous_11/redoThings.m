function redoThings

clear
clc

% ======================= Start the counter ===============================

tic
% counter_init=cputime;
counter_init=clock;

% ===================== User defined variables ============================

nex=2;
ney=2;
nez=3;

dto=2.0e-6;
dt=2.0e-6;

obs_t=2.0e-1;

diff=1800;

T=0.6;

ci = 0.7; %Overall initial concentrations
ci_fluc = 0.01; %Initial fluctuation

nr_tol=1.0e-6;
nr_max_iteration=20;

n1=1;
n2=10;

entropy=1;
T_theta=1;

xlen=1;
ylen=1;
zlen=1;

fps=10;

dt_down=0.5;
dt_up=1.2;
dt_ideal=dt;
bypass_tol=5;

time_out=3600*8;
% time_out=600;

max_frame=10000;

% === Specify constants which will be used throughout =====================

w=[0.27778 0.4444 0.27778];
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

nx=nex+1;ny=ney+1;nz=nez+1;
ne=nex*ney*nez;
n=nx*ny*nz;
neight=8*n;

ll=2;li=8;lu=2+8*ney;
bl=3;bi=8*ny;bu=3+8*nex*ny;
tl=3+8*ney;ti=8*ny;tu=3+8*ney+8*nex*ny;
rl=2+8*nex*ny;ri=8;ru=2+8*ney+8*nex*ny;
zi=nx*ny*8;zu=zi*(nz-1);bny=nx*ny*8;

weights=set_up_rec_mesh_3D(gp);

% [irows,icols,isj,col_placement,nnz_]=compact_rec_3D(n,nx,ny,nz);
% hello=1:1:nnz_
% iSparse=sparse(irows,icols,hello)

% iSparse=zeros(neight,neight);
% for i=1:1:nnz_
%     iSparse(irows(i),icols(i))=i;
% end

% === initialize time =====================================================

time=0.0;

% === Evaluate value of chi ===============================================

chi=0.5-entropy*(1-T_theta/T);

% === Determine the coordinate at each node ===============================

[x,y,z,x_coord,y_coord,z_coord]=...
    rectangular_meshing_3D(nex,ney,nez,nx,ny,nz,n,xlen,ylen,zlen);

% === Randomize and apply initial conditions ==============================

co=generate_co(ci,ci_fluc,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny,neight);
ci_ave=mean2(extract_abs_result_rec_3D(co,nx,ny,nz));

% === Take the very first time step =======================================

c=co;

% === Set up for plotting =================================================

% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% fig=figure(1);
% set(fig,'Position',[0 0 1280 1080])
% minz=0.0;maxz=1.0;
% c_obs=[ci-ci_fluc ci_ave ci+ci_fluc];
% face_colors={'blue','yellow','red'};
% view_angle=3;
% transparency=0.9;
% edge_color='none';
% fri=0;
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID

% === Plot the concentration in second step, which is qual to first one ===
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% fri=fri+1;
% clf(fig)
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% visual(fri)=visual3D(c,c_obs,x_coord,y_coord,z_coord,...
%     face_colors,edge_color,nx,ny,nz,view_angle,transparency,time);

% visual(fri)=visual_surf_3D(c,nx,ny,nz,x_coord,y_coord,z_coord);

% visual(fri)=visual_rec_3D(c,ci_ave,nx,ny,nz,x_coord,y_coord,z_coord);

% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% % HYBRID
% visual(fri)=visual_all_3D(c,c_obs,x_coord,y_coord,z_coord,...
%     face_colors,edge_color,nx,ny,nz,view_angle,transparency,time,minz,maxz);
% HYBRID
% HYBRID
% HYBRID
% HYBRID
% HYBRID
% HYBRID
% HYBRID
% error('Just stop')

% === Prepare to enter main loop ==========================================

bypass=0;
bypass_count=0;

% === THESE THINGS HERE CAN BE SUMMED UP INTO SETTING UP PROCESS LATER! ===


% error('Just stop')
%==========================================================================

%====================================== Test
[gbfs1,gbfs2]=sj_mapper_3d(nx,ny,nz);
nnz_=nnz_sj_3d(nx,ny,nz);
weights222=weight_adjuster_3d(generateWeights_3d(),1/nex,1/ney,1/nez);
% conc=get_conc_3d(ne,c,weights222,nex*ney,nex,nx,ny);
% test222=compute_sj_3d(gbfs1,gbfs2,nnz_,conc,dt,weights222,nex,ney,nx,ny,nz,n1,n2,1/(nex*ney*nez));

test101=generateWeights_3d();

[1/nex,1/ney,1/nez]

size(weights)



[test101(1,1,3,1,2,1) weights222(1,1,3,1,2,1) weights(1,3,4)]

% get_conc_3d(ne,c,weights,nex*ney,nex,nx,ny)

disp('hello count the stars')
% return

test_phi=zeros(64,9);
itest=0;
%========================================

% for rahmat=1:1:10
while 1
    
    cp=c+dt*(c-co)./dto;
    coo=co;
    co=c;
    c=cp;
    c=bc3D(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny);
    
    err=nr_tol+1;
    jk=0;
    
    
    %======================================== Test
    disp('Computing new method')
    conc=get_conc_3d(ne,c,weights222,nex*ney,nex,nx,ny);
    test_sj=compute_sj_3d(gbfs1,gbfs2,nnz_,conc,dt,weights222,nex,ney,nx,ny,nz,n1,n2,1/(nex*ney*nez));
    disp('End new method')
    %========================================
    
    %for saeedi=1:10
    while err > nr_tol
        sf=zeros(1,neight);
        sj=zeros(neight,neight);
%         sj_compact=zeros(1,nnz_);
        fprintf('Filling in sf sj ')
        
        
        for e=1:1:ne
%             e
            [gbasis,dx,dy,dz]=echaracterize(e,nex,ney,nez,nx,ny,nz,x,y,z);
            iweight=0;
            for w1=1:1:3
                for w2=1:1:3
                    for w3=1:1:3
                        iweight=iweight+1;
                        [phi,phix,phiy,phiz,phixx,phiyy,phizz]=...
                            extract_weights(iweight,weights);
                        [phi,phix,phiy,phiz,phixx,phiyy,phizz]=...
                            rec_basis_transform_3D(dx,dy,dz,phi,...
                            phix,phiy,phiz,phixx,phiyy,phizz);
                        [con,cono,conx,cony,conz,conxx,conyy,conzz]=...
                            eval_local3D(c,co,phi,phix,phiy,phiz,...
                            phixx,phiyy,phizz,gbasis);
                        
                        
                        %=========================================================== Test
                        if e==3 %&& w1==1 && w2==2 && w3==3
                            
%                             conc(e,1,w3,w2,w1)
%                             conc(e,1,w1,w2,w3)
%                             con
%                             error('just stop here')
%                             size(phi)
%                             test_phi(w3,w2,w1)=phi;
                            itest=itest+1;
                            test_phi(:,itest)=phi';
                            
                        end
                        %===========================================================
                        
                        
                        [sf,sj]=eval_sfsj_3D(sf,sj,w(w1),w(w2),w(w3),...
                            dx,dy,dz,n1,n2,con,cono,conx,cony,conz,...
                            conxx,conyy,conzz,phi,phixx,phiyy,phizz,...
                            gbasis,dt,diff,T,chi);
%                         [sf,sj_compact]=eval_compactly_sfsj_3D(sf,sj_compact,isj,...
%                             iSparse,col_placement,w(w1),w(w2),w(w3),...
%                             dx,dy,dz,n1,n2,con,cono,conx,cony,conz,...
%                             conxx,conyy,conzz,phi,phixx,phiyy,phizz,...
%                             gbasis,dt,diff,T,chi);
                    end
                end
            end
        end
        %===========================================
% %         [irows,icols]=compact_rec_3D(n,nx,ny,nz);
% %         [iirows,iicols]=find(sj);
% %         length(irows)
% %         length(iirows)
% %         length(sj_compact)
% %         nnz_
% %         sj_compact;
% %         nnz(sj_compact)
%         dxyz=dx*dy*dz;
% 
%         weights222=weight_adjuster_3d(generateWeights_3d(),1/nex,1/ney,1/nez);
%         conc_=get_conc_3d(ne,c,weights222,nex*ney,nex,nx,ny);
%         conco_=conc_;
%         sf_=compute_sf_3d(neight,conc_,conco_,nx,ny,nz,nex,ney,weights222,dt,n1,n2,dxyz);
%         
%         [sf' sf_']
% 
%         disp('end test here------------------')
%         return
% 
%         sj;
%         
%         con;
%         
%         [gbfs1,gbfs2]=sj_mapper_3d(nx,ny,nz);
%         nnz_=nnz_sj_3d(nx,ny,nz);
%         weightsOriginal=weights;
%         weights=generateWeights_3d();
%         
%         size(weightsOriginal);
%         size(weights);
%         
%         sum(sum(sum(weightsOriginal)));
%         
%         sum(sum(sum(sum(sum(sum(weights))))));
%         
% %         return
%         
%         
%         weights=weight_adjuster_3d(generateWeights_3d(),1/nex,1/ney,1/nez);
%         
%         
%         
%         
%         
% %         return
%         
%         conc=get_conc_3d(ne,c,weights,nex*ney,nex,nx,ny);
%         
% %         test;
% %         
% %         size(test)
% %         
% %         nnz(test)
%         
%         test=compute_sj_3d(gbfs1,gbfs2,nnz_,conc,dt,weights,nex,ney,nx,ny,nz,n1,n2,1/(nex*ney*nez));
%         
%         sj;
%         
%         test101=zeros(1,nnz_);
%         
% %         sj_test=zeros(neight,neight);
%         for i=1:1:nnz_
%             test101(i)=sj(gbfs1(i),gbfs2(i));
%             if test101(i)==0
%                 error('just sttestesdf')
%             end
%         end
%         
% %         test=test*dx*dy*dz;
% %         
%         diff=test-test101;
%         
%         disp('Test    Correct   Difference')
%         
%         A=[test222' test101' diff'];
%         
%         disp(A(1:10,:))
%         
%         max(diff);
%         
%         
%         
%         
%         return
%         %-------
%         
%         taa=get_gbfs_of_element(1,nex*ney,nex,nx,ny);
%         takeit=zeros(1,64);
%         for honey=1:1:64
%             takeit(honey)=c(taa(honey));
%         end
%         
%         takeit'
%         
%         A(1,1)
%         
%         dx
%         1/nex
%         
%         error('n\Good')
% === Apply boundary conditions to sf and sj ==============================

        if time~=0
            nnz(sj)
            error('Just stop')
        end

        fprintf('Applying BC ')
        [sf,sj]=bsfsfsj3D(sf,sj,...
            ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,bny,nx,ny,nz,neight);
%         error('Just stop')
% === Create compact matrix ===============================================
        fprintf('Compressing ')
        toc
        [rowi,coli,val]=find(sj);
        toc
        fprintf('Sparsing ')
        sj=sparse(rowi,coli,val);
        toc
        
% === Carry matrix Division ===============================================
        
        fprintf('Dividing ')
        c_=sj\sf';
        toc
        
% === Update the solution =================================================
        fprintf('Updating the solution ')
        c = c + c_';
        
% === Evaluate the error ==================================================
        fprintf('Evaluate the error ')
        err = sqrt(sum(c_.^2.0))
        
% === Assert reasonable convergence =======================================

        jk = jk + 1;
        if (jk>1 && err>=erro) || jk>=nr_max_iteration
            fprintf('Newton-Raphson is not converging ')
            jk=-1;
            break
        else
            erro=err;
        end
    end
    
    fprintf('Newton-Raphson cycle finishes ')
    
% === Resolve divergence ==================================================

    if jk==-1
        bypass=bypass+1;
        if bypass>bypass_tol
            fprintf('Too much bypassing!')
            break
        else
            c=co;
            co=coo;
            dt=dt*dt_down
            bypass_count=bypass_count+1;
            continue
        end
    elseif bypass~=0
        bypass=0;
    end
    
% === Increment time step =================================================
    
    time=time+dt % Time for last loop, not next!
    
% === Visualize the result ================================================
    
    fri=fri+1;
    clf(fig)
%     visual(fri)=visual3D(c,c_obs,x_coord,y_coord,z_coord,...
%     face_colors,edge_color,nx,ny,nz,view_angle,transparency,time);

%     visual(fri)=visual_rec_3D(c,ci_ave,nx,ny,nz,x_coord,y_coord,z_coord);
    visual(fri)=visual_all_3D(c,c_obs,x_coord,y_coord,z_coord,...
        face_colors,edge_color,nx,ny,nz,view_angle,transparency,time,...
        minz,maxz);
    if fri==max_frame
        fprintf('Max frame reached!')
        break
    end
    
% === Time out check point ================================================

    if time_out<etime(clock,counter_init)
        fprintf('Time out!')
        break
    end
    
% === Prepare for next time loop ==========================================

    if time>=obs_t
        fprintf('Observation time goal reached!')
        break
    end
    
    dtoo=dto;
    dto=dt;
    
    if jk>=10
        dt=dt*dt_down;
    elseif jk<=3 || dt<dt_ideal
        dt=dt*dt_up;
    end
    
    if time>obs_t
        time=obs_t;
        dt=obs_t-(time-dt);
    end
    
end

% === Export the movie ====================================================

video=VideoWriter('Nonlinear Cahn Hilliard 3D.avi', 'Uncompressed AVI');
video.FrameRate=fps;
open(video)
writeVideo(video,visual);
close(video)

% === Conclusion ==========================================================

toc
bypass_count
movie(gcf,visual,1)

end