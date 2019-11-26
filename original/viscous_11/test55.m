function test55

clear
clc

nex=5;
ney=4;
nez=3;

%--- Set Up --------------------------------

nx=nex+1;
ny=ney+1;
nz=nez+1;

n=nx*ny*nz;

nxny=nx*ny;
nexney=nex*ney;

ne=nex*ney*nez;

%---Unique for sphere---
nx_=nx-2;
ny_=ny-2;
nz_=nz-2;
n_nodesurf=2*(nx_*ny_+nx_*nz_+ny_*nz_)+4*(nx_+ny_+nz_)+8;

n_nodessurf_x=ny*nz;
n_nodessurf_y=nx_*nz;
n_nodessurf_z=nx_*ny_;

n_nSurf=zeros(1,6);
n_nSurf(1)=n_nodessurf_x;
n_nSurf(2)=2*n_nodessurf_x;
n_nSurf(3)=n_nSurf(2)+n_nodessurf_y;
n_nSurf(4)=n_nSurf(3)+n_nodessurf_y;
n_nSurf(5)=n_nSurf(4)+n_nodessurf_z;
n_nSurf(6)=n_nSurf(5)+n_nodessurf_z;

if n_nSurf(6)~=n_nodesurf
    error('dafafafdagasd')
end

n_eSurf=zeros(1,6);

n_eSurf_x=ney*nez;
n_eSurf_y=nex*nez;
n_eSurf_z=nex*ney;

n_eSurf(1)=n_eSurf_x;
n_eSurf(2)=n_eSurf_x*2;
n_eSurf(3)=n_eSurf(2)+n_eSurf_y;
n_eSurf(4)=n_eSurf(3)+n_eSurf_y;
n_eSurf(5)=n_eSurf(4)+n_eSurf_z;
n_eSurf(6)=n_eSurf(5)+n_eSurf_z;

n_eSurf;


neight=8*n;

neightTotal=neight+8*n_nSurf(6);

nTotal=n+n_nSurf(6);

neightTotal=8*nTotal;

%-----------------------------

%--------------------------------------

% % [X,Y,Z]=generate_sphere_mesh_with_ratio(nx,ny,nz,n);
[X,Y,Z]=generate_sphere_mesh_with_extrusion_3d(nx,ny,nz,n,n_eSurf,n_nSurf,ne);
% 
% 
% 
% e=136
% nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_)
% 
% 
% node=44;
% a=[X(node) Y(node) Z(node)]
% b=sqrt(sum(a.^2))

%--- Plotting ------------------------------------


% plot3(X,Y,Z,'.')
% grid on
% xlabel('x');ylabel('y');zlabel('z');

%--------------------------------------------------

%---Get nodes of elements----------------------------

% connectivity=...
%     [1 2
%     3 4
%     1 3
%     2 4
%     5 6
%     7 8
%     5 7
%     6 8
%     1 5
%     2 6
%     3 7
%     4 8];
% hold on
% for e=1:1:ne
%     nodes=get_nodes_of_elements_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%     for i=1:1:12
%         xLine=[X(nodes(connectivity(i,1))) X(nodes(connectivity(i,2)))];
%         yLine=[Y(nodes(connectivity(i,1))) Y(nodes(connectivity(i,2)))];
%         zLine=[Z(nodes(connectivity(i,1))) Z(nodes(connectivity(i,2)))];
%         plot3(xLine,yLine,zLine,'g')
%     end
% end
% 
% for e=ne+1:1:n_eSurf(6)+ne
% % for e=n_eSurf(5)+ne+1:1:n_eSurf(6)+ne
%     nodes=get_nodes_of_elements_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%     for i=1:1:12
%         xLine=[X(nodes(connectivity(i,1))) X(nodes(connectivity(i,2)))];
%         yLine=[Y(nodes(connectivity(i,1))) Y(nodes(connectivity(i,2)))];
%         zLine=[Z(nodes(connectivity(i,1))) Z(nodes(connectivity(i,2)))];
%         plot3(xLine,yLine,zLine,'r')
%     end
% end
% hold off
%--------------------------------------------------

%---Extrusion------------------------------------
% hold on
% for node=1:1:n
%     [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
%     if xth==1 || xth==nx || yth==1 || yth==ny || zth==1 || zth==nz
% 
%         ext=get_extrusion_sph(node,n,nx,ny,nz,n_nSurf,nx_);
%         x=[X(node) X(ext)];
%         y=[Y(node) Y(ext)];
%         z=[Z(node) Z(ext)];
%         plot3(x,y,z,'b')
%     end
% end
% hold off
%--------------------------------------------------

%---Adjacency-----------------------------------------------

% node=192
% 
% [nodes_adjacent,extruding]=get_node_adjacent_sph(node,n,nx,ny,nz,nxny,n_nSurf,nx_)
% 
% len=length(extruding);
% hold on
% plot3(X(node),Y(node),Z(node),'r*')
% for i=1:1:len
%     node=extruding(i);
%     plot3(X(node),Y(node),Z(node),'ro')
% end
% hold off

%--------------------------------------------------

%---Extrusion inverse-----------------------------------------------
% hold on
% for node=n+1:1:n_nSurf(6)+n
%     [xth,yth,zth]=get_root_of_extrusion_sph(node,n_nSurf,n,nx,ny,nz,nx_);
%     root=get_node_3d(xth,yth,zth,nx,ny);
%     x=[X(node) X(root)];
%     y=[Y(node) Y(root)];
%     z=[Z(node) Z(root)];
%     plot3(x,y,z,'r')
% end
% hold off
%--------------------------------------------------

%---Extruding element-----------------------------------------------

% e=57
% get_extruding_element_sph(e,nex,ney,nez,nexney,n_eSurf)

%--------------------------------------------------

%---Mirror-----------------------------------------------

% hey=mirror(1,1)

%--------------------------------------------------

%---Analyzing interactions-----------------------------------------------

% gbf2=1697;
% gbf1=1704;
% 
% [elements,lbfs1,lbfs2,mx1,mx2,mx3,mx4,mx5,mx6]=...
%     analyze_interaction_sph(gbf1,gbf2,ne,nex,ney,nez,n,nx,ny,nz,neight,...
%     nexney,n_eSurf,n_nSurf,nx_);
% 
% 
% 
% [elements' lbfs1' lbfs2']
% 
% mx1
% mx2
% mx3
% mx4
% mx5
% mx6
% 
% %Come up with algorithm to verify!!!

%--------------------------------------------------

%---nnz of new sj-----------------------------------------------
% 
% disp('yo')
% 
% nnz_=nnz_sj_sph(nx,ny,nz,n_nSurf)

%--------------------------------------------------

%---sj mapper-----------------------------------------------

% [gbfs1,gbfs2]=sj_mapper_sph(n,nx,ny,nz,nnz_,nxny,n_nSurf,nx_,nTotal);
% 
% nnz(gbfs1)
% nnz(gbfs2)

%--------------------------------------------------

%---gbfs of elements-----------------------------------------------

% gbfs_=get_gbfs_of_element_sph(5,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%--------------------------------------------------

%---Verification------------------------------------
% tic
% parfor i=1:1:nnz_
%     [elements,lbfs1,lbfs2,mx1,mx2,mx3,mx4,mx5,mx6]=...
%         analyze_interaction_sph(gbfs1(i),gbfs2(i),ne,nex,ney,nez,n,nx,ny,nz,neight,...
%         nexney,n_eSurf,n_nSurf,nx_);
%     
%     %--------------------------------------------------
%     if ~isnan(elements)
%         for ie=1:1:8
%             if elements(ie)~=0
%                 e=elements(ie);
%                 lbf1=lbfs1(ie);
%                 lbf2=lbfs2(ie);
%     %             nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
% 
%                 gbfs_=get_gbfs_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
% 
% 
%                 compare1=gbfs_(lbf1);
%                 compare2=gbfs_(lbf2);
% 
%                 if compare1~=gbfs1(i) || compare2~=gbfs2(i)
%                     error('bad')
%                 end
%             end
%         end
%     end
%     %---------------------------------------------
%     
%     mx=mx6;
%     if ~isnan(mx)
%         test=mx(:,1);
%         if sum(test)==0
%             error('so bad')
%         end
%         for ie=1:1:4
%             if mx(ie,1)~=0
%                 e=mx(ie,1);
%                 lbf1=mx(ie,2);
%                 lbf2=mx(ie,3);
%     %             nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%                 gbfs_=get_gbfs_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%                 compare1=gbfs_(lbf1);
%                 compare2=gbfs_(lbf2);
% 
%                 if compare1~=gbfs1(i) || compare2~=gbfs2(i)
%                     error('bad')
%                 end
%             end
%         end
%     end
%     
%     
%     %---------------------------------------------
%     
%     
%     
% end
% toc
%----------------------------------------------------

%---Adjacent elements------------------------------------
% 
% get_extruding_element_sph(57,nex,ney,nez,nexney,n_eSurf)
% 
% 
% get_esurf_minus_x(ney,2,1)

% node=202
% 
% for node=1:1:nTotal
%     [elements,mx1,mx2,mx3]=get_adjacent_elements_sph(node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_);
% end
% %VERIFY THIS LATER!!!!!!!!!!!!!!!!!
% 
% for node=1:1:n
%     [elements,mx1,mx2,mx3]=get_adjacent_elements_sph(node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_);
%     [xth,yth,zth]=get_xyzth_3d(node,nx,ny);
% 
%     if xth==1 || xth==nx || yth==1 || yth==ny || zth==1 || zth==nz
%         sum1=0;
%         if ~isnan(mx1)
%             sum1=sum1+nnz(mx1);
%         end
%         if ~isnan(mx2)
%             sum1=sum1+nnz(mx2);
%         end
%         if ~isnan(mx3)
%             sum1=sum1+nnz(mx3);
%         end
%         
%         if (xth==1 || xth==nx) && (yth==1 || yth==ny) && (zth==1 || zth==nz)
%             
%             if sum1~=3
%                 error('bad')
%             end
%         else
% %             warning('dafadf')
%             if sum1~=4
%                 error('bad')
%             end
%         end
%     end
% end
% 
% 
% for node=1:1:nTotal
%     [elements,mx1,mx2,mx3]=get_adjacent_elements_sph(node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_);
% end
% %VERIFY THIS LATER!!!!!!!!!!!!!!!!!
% 
% for node=1+n:1:nTotal
%     [elements,mx1,mx2,mx3]=get_adjacent_elements_sph(node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_);
%     [xth,yth,zth]=get_root_of_extrusion_sph(node,n_nSurf,n,nx,ny,nz,nx_);
% 
%     if xth==1 || xth==nx || yth==1 || yth==ny || zth==1 || zth==nz
%         sum1=0;
%         if ~isnan(mx1)
%             sum1=sum1+nnz(mx1);
%         end
%         if ~isnan(mx2)
%             sum1=sum1+nnz(mx2);
%         end
%         if ~isnan(mx3)
%             sum1=sum1+nnz(mx3);
%         end
%         if (xth==1 || xth==nx) && (yth==1 || yth==ny) && (zth==1 || zth==nz)
%             
%             if sum1~=3
%                 error('bad')
%             end
%         else
% %             warning('dafadf')
%             if sum1~=4
%                 error('bad')
%             end
%         end
%     end
% end

%----------------------------------------------------

%==============================
%==============================
%==============================
%==============================
%==============================
%==============================

% %---Get nodes of elements----------------------------
% 
% connectivity=...
%     [1 2
%     3 4
%     1 3
%     2 4
%     5 6
%     7 8
%     5 7
%     6 8
%     1 5
%     2 6
%     3 7
%     4 8];
% hold on
% for e=1:1:ne
%     nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%     for i=1:1:12
%         xLine=[X(nodes(connectivity(i,1))) X(nodes(connectivity(i,2)))];
%         yLine=[Y(nodes(connectivity(i,1))) Y(nodes(connectivity(i,2)))];
%         zLine=[Z(nodes(connectivity(i,1))) Z(nodes(connectivity(i,2)))];
%         plot3(xLine,yLine,zLine,'g')
%     end
% end
% 
% for e=ne+1:1:n_eSurf(6)+ne
% % for e=n_eSurf(5)+ne+1:1:n_eSurf(6)+ne
%     nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
%     for i=1:1:12
%         xLine=[X(nodes(connectivity(i,1))) X(nodes(connectivity(i,2)))];
%         yLine=[Y(nodes(connectivity(i,1))) Y(nodes(connectivity(i,2)))];
%         zLine=[Z(nodes(connectivity(i,1))) Z(nodes(connectivity(i,2)))];
%         plot3(xLine,yLine,zLine,'r')
%     end
% end
% hold off
% grid on
%--------------------------------------------------

%===SEPTEMBER VER PART================================

nnz_=nnz_sj_sph(nx,ny,nz,n_nSurf);
[gbfs1,gbfs2]=sj_mapper_sph(n,nx,ny,nz,nnz_,nxny,n_nSurf,nx_,nTotal);

counter=zeros(1,neightTotal);

for sjth=nnz_
    
%     if mod(sjth,100)==0
%         [sjth nnz_]
%     end
    
    gbf1=gbfs1(sjth);
    gbf2=gbfs2(sjth);
%     
%     gbf1=7
%     gbf2=1537
%     
% %     if gbf1==7 && gbf2==1537
% %         error('ok')
% %     end
%     
%     [elements,lbfs1,lbfs2,mx1,mx2,mx3,mx4,mx5,mx6]=...
%         analyze_interaction_sph(gbf1,gbf2,ne,nex,ney,nez,n,nx,ny,nz,neight,...
%         nexney,n_eSurf,n_nSurf,nx_)
%     counter(gbf1)=counter(gbf1)+1;
%     if gbf1==

%     for ijk=1:1:nnz_
%         if gbfs1(ijk)==gbf1 && gbfs2(ijk)==gbf2 && sjth~=ijk
%             error('overlap')
%         end
%     end

end

% nodes=get_nodes_of_element_sph(115,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
% nodes'
% 
% 193*8

% counter=counter/8;
% 
% counter';
% max(gbfs1);
% counter(1280)

for node=1:1:nTotal
    [elements,eExtruding_x,eExtruding_y,eExtruding_z,x_dir,y_dir,z_dir]=...
        get_adjacent_elements_sph(node,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_);
%     if ~isnan(elements)
%         size(elements);
%     end
    mx=elements;
    if ~isnan(mx)
        for ie=1:1:4
            if mx(ie)~=0
                e=mx(ie);
                nodes=get_nodes_of_element_sph(e,ne,nex,ney,nexney,nx,ny,nz,n,nxny,n_eSurf,n_nSurf,nx_);
                if ~ismember(node,nodes)
                    error('love')
                end
            end
        end
    end
end


[elements,eExtruding_x,eExtruding_y,eExtruding_z,x_dir,y_dir,z_dir]=...
    get_adjacent_elements_sph(9,n,nx,ny,nz,ne,nex,ney,n_eSurf,n_nSurf,nTotal,nx_)


end



function sol=get_esurf_minus_x(ney,yth,zth)
sol=(zth-1)*ney+yth;
end

function sol=get_esurf_plus_x(n_eSurf,ney,yth,zth)
sol=n_eSurf(1)+(zth-1)*ney+yth;
end

function sol=get_esurf_minus_y(n_eSurf,nex,xth,zth)
sol=n_eSurf(2)+(zth-1)*nex+xth;
end

function sol=get_esurf_plus_y(n_eSurf,nex,xth,zth)
sol=n_eSurf(3)+(zth-1)*nex+xth;
end

function sol=get_esurf_minus_z(n_eSurf,nex,xth,yth)
sol=n_eSurf(4)+(yth-1)*nex+xth;
end

function sol=get_esurf_plus_z(n_eSurf,nex,xth,yth)
sol=n_eSurf(5)+(yth-1)*nex+xth;
end