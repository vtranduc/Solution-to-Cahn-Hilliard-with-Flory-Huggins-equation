function [sf,sj]=bsfsfsj3D(sf,sj,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,zi,zu,nxy,nx,ny,nz,neight)

nxny=nx*ny;

for i=ll:li:lu
    for j=i:zi:i+zu
        sf(j)=0.0;
%         for k=1:1:neight
%             sj(j,k)=0.0;
%         end
%         sj(j,j)=1.0;
%         node=(j-2)/8+1;
%         interactions=interactions_rec_3D(node,ny,nz,nxny);
%         for k=1:1:length(interactions)
%             pre_index=interactions(k)*8-8;
%             for l=1:1:8
%                 sj(j,pre_index+l)=0.0;
%             end
%         end
        sj=specific_row_reset_rec_3D(sj,j,ny,nz,nxny);
        sj(j,j)=1.0;
    end
end
for i=bl:bi:bu
    for j=i:zi:i+zu
        sf(j)=0.0;
%         for k=1:1:neight
%             sj(j,k)=0.0;
%         end
        sj=specific_row_reset_rec_3D(sj,j,ny,nz,nxny);
        sj(j,j)=1.0;
    end
end
for i=rl:ri:ru
    for j=i:zi:i+zu
        sf(j)=0.0;
%         for k=1:1:neight
%             sj(j,k)=0.0;
%         end
        sj=specific_row_reset_rec_3D(sj,j,ny,nz,nxny);
        sj(j,j)=1.0;
    end
end
for i=tl:ti:tu
    for j=i:zi:i+zu
        sf(j)=0.0;
%         for k=1:1:neight
%             sj(j,k)=0.0;
%         end
        sj=specific_row_reset_rec_3D(sj,j,ny,nz,nxny);
        sj(j,j)=1.0;
    end
end

for i=4:8:nxy-4
    sf(i)=0.0;
    sf(i+zu)=0.0;
%     for k=1:1:neight
%         sj(i,k)=0.0;
%         sj(i+zu,k)=0.0;
%     end
    sj=specific_row_reset_rec_3D(sj,i,ny,nz,nxny);
    sj=specific_row_reset_rec_3D(sj,i+zu,ny,nz,nxny);
%     for kkk=1:1:neight
%         if sj(i,kkk)~=0 || sj(i+zu,kkk)~=0
%             error('just raise')
%         end
%     end
    sj(i,i)=1.0;
    sj(i+zu,i+zu)=1.0;
    
end
% sj(i+zu,:)
end