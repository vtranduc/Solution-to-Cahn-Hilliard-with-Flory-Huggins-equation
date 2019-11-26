function build_3D_inProcess

clear
clc
tic
% nex=5;
% ney=4;
% nez=3;

nex=5;
ney=4;
nez=3;

nx=nex+1;
ny=ney+1;
nz=nez+1;

disp('nnz is')

nnz_=nnz_sj_3d(nx,ny,nz)

node=75;

%====Verified!!!=====
[xth,yth,zth]=get_xyzth_3d(node,nx,ny);
get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);
%====================

% node=get_node_3d(2,2,1,nx,ny);
% 
% gbs_nodal=get_gbs_3d(66);
% 
% test=get_node_adjacent_3d(1,nx,ny,nz);
% 
% test';
% 
% [gbs1,gbs2]=sj_mapper_3d(nx,ny,nz);
% 
% wah=[gbs1' gbs2'];
% 
% nnz(gbs2)
% 
% size(gbs2);
% 
% size(gbs1)
% 
% nnz_
% 
% test101=get_elements_3d(3,4,2,nex,ney,nx,ny,nz);
% 
% test101';
% 
% [node,type]=analyze_gbs_3d(601)

% analyze_interaction_3d(530,290,nx,ny,nz)

% dlmwrite('show_me.txt',wah)

[gbfs1,gbfs2]=sj_mapper_3d(nx,ny,nz);

nnz_

length(gbfs1)

length(gbfs2)

% solution=get_elements_3d(2,2,2,nex,ney,nx,ny,nz,...
%     3)
% return
% 
% for i=1:1:nnz_
%     analyze_interaction_3d(gbfs1(i),gbfs2(i),nex,ney,nx,ny,nz);
% end

[elements,lbfs1,lbfs2]=analyze_interaction_3d(353,311,nex,ney,nx,ny,nz);

[elements',lbfs1',lbfs2'];

nexney=nex*ney;

solution=get_gbfs_of_element(1,nexney,nex,nx,ny);

solution'

toc
return
%=================================
%=================================
%=================================
%=================================
%=================================
%=================================

[xth,yth,zth,nodeth,relationshipth]=get_xyzth_sj_3d(133120,nx,ny,nz,nnz_) %uncertain


tester1=zeros(nnz_,4);
tester100=zeros(nnz_,5);
for i=1:1:nnz_
    [xth,yth,zth,nodeth,relationth]=get_xyzth_sj_3d(i,nx,ny,nz,nnz_); %uncertain
%     tester1(i,:)=get_xyzth_sj_3d(i,nx,ny,nz,nnz_);
    tester1(i,:)=[xth,yth,zth,relationth];
    tester100(i,:)=[xth,yth,zth,nodeth,relationth];
end
tester1(1:200,:)

tester=0;
xp=1;
yp=1;
zp=1;
tester2=0;
for i=1:1:nnz_
    if tester1(i,4)~=tester+1
        %check
        if tester1(i,4)~=1
%             tester
%             tester1(i-2,:)
%             tester1(i-1,:)
%             tester1(i,:)
%             warning('bad')
        end
        if tester1(i,1)-xp==1 || tester1(i,2)-yp==1 || tester1(i,3)-zp==1
%             disp('change!')
            xp=tester1(i,1);
            yp=tester1(i,2);
            zp=tester1(i,3);
            tester2=tester2+tester1(i-1,4);
            tester1(i-1,4);
            tester=1;
        else
%             warning('Bad')
        end
    else
        tester=tester+1;
    end
end
tester2=tester2+tester1(nnz_,4);

tester2
nnz_

tester1;


tester100

% dlmwrite('just_output.txt',tester100)

end