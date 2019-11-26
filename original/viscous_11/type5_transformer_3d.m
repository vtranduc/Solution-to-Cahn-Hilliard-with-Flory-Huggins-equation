function sol=type5_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol)

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation=0;
                for ibeta=1:1:3
                    summation=summation+w(ibeta)...
                        *local_basis_3d(gps(ix),gps_plus(iy,ibeta),gps(iz),...
                        orientation_list(iorientation,:),[1 1 0],[1 1 0]);
                end
                sol(iorientation,5,2,ix,iy,iz)=summation*gps(iy)...
                    *map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation=0;
                for ialpha=1:1:3
                    summation=summation+w(ialpha)...
                        *local_basis_3d(gps_plus(ix,ialpha),gps(iy),gps(iz),...
                        orientation_list(iorientation,:),[1 1 0],[1 1 0]);
                end
                sol(iorientation,5,3,ix,iy,iz)=summation*gps(ix)...
                    *map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 0 0]);
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation1=0;
                summation2=0;
                summation3=0;
                for ialpha=1:1:3
                    for ibeta=1:1:3
                        term1=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eX,[1 0 0]);
                        term2=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eY,[0 1 0]);
                        term3=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eX,[0 1 0]);
                        term4=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eY,[1 0 0]);
                        determinant=term1*term2-term3*term4;
                        fz=get_fg_dervs_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),...
                            orientation_list(iorientation,:),[1 1 0],eX,eY,eZ,3);
                        summation1=summation1+w(ialpha)*w(ibeta)*determinant*fz;
                        fzz=get_fg_dervs_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),...
                            orientation_list(iorientation,:),[1 1 0],eX,eY,eZ,6);
                        summation2=summation2+w(ialpha)*w(ibeta)*determinant*fzz;
                        summation3=summation3+w(ialpha)*w(ibeta)*determinant...
                            *local_basis_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),...
                            orientation_list(iorientation,:),[1 1 0],[1 1 0]);
                    end
                end
                sol(iorientation,5,4,ix,iy,iz)=summation1*gps(ix)*gps(iy);
                sol(iorientation,5,7,ix,iy,iz)=summation2*gps(ix)*gps(iy);
                sol(iorientation,5,1,ix,iy,iz)=summation3*gps(ix)*gps(iy);
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation=0;
                for ibeta=1:1:3
                    fx=get_fg_dervs_3d(gps(ix),gps_plus(iy,ibeta),gps(iz),...
                        orientation_list(iorientation,:),[1 1 0],eX,eY,eZ,1);
                    summation=summation+w(ibeta)*fx;
                end
                sol(iorientation,5,5,ix,iy,iz)=summation*gps(iy)...
                    *map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation=0;
                for ialpha=1:1:3
                    fy=get_fg_dervs_3d(gps_plus(ix,ialpha),gps(iy),gps(iz),...
                        orientation_list(iorientation,:),[1 1 0],eX,eY,eZ,2);
                    summation=summation+w(ialpha)*fy;
                end
                sol(iorientation,5,6,ix,iy,iz)=summation*gps(ix)...
                    *map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 0 0]);
            end
        end
    end
end

end