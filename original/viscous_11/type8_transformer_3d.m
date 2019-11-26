function sol=type8_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol)

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %----------------
                summation=0;
                J=zeros(3,3);
                for ialpha=1:1:3
                    for ibeta=1:1:3
                        for igamma=1:1:3
                            J(1,1)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eX,[1 0 0]);
                            J(1,2)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eX,[0 1 0]);
                            J(1,3)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eX,[0 0 1]);
                            J(2,1)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[1 0 0]);
                            J(2,2)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 1 0]);
                            J(2,3)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 0 1]);
                            J(3,1)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[1 0 0]);
                            J(3,2)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 1 0]);
                            J(3,3)=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 0 1]);
                            summation=summation+w(ialpha)*w(ibeta)*w(igamma)*det(J)...
                                *local_basis_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps_plus(iz,igamma),...
                                orientation_list(iorientation,:),[1 1 1],[1 1 1]);
                        end
                    end
                end
                sol(iorientation,8,1,ix,iy,iz)=summation*gps(ix)*gps(iy)*gps(iz);
                %-----------------------
            end
        end
    end
end

%=====================================
%=====================================

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %------------------------------------
                summation=0;
                for ibeta=1:1:3
                    for igamma=1:1:3

                        term1=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 1 0]);
                        term2=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 0 1]);
                        term3=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 0 1]);
                        term4=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 1 0]);
                        determinant=term1*term2-term3*term4;
                        summation=summation+w(ibeta)*w(igamma)*determinant...
                            *local_basis_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),orientation_list(iorientation,:),[1 1 1],[1 1 1]);
                    end
                end
                sol(iorientation,8,2,ix,iy,iz)=summation*gps(iy)*gps(iz);
                %-------------------------------------
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %------------------------------------
                summation=0;
                for ialpha=1:1:3
                    for igamma=1:1:3
                        term1=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eX,[1 0 0]);
                        term2=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eZ,[0 0 1]);
                        term3=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eX,[0 0 1]);
                        term4=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eZ,[1 0 0]);
                        determinant=term1*term2-term3*term4;
                        summation=summation+w(ialpha)*w(igamma)*determinant...
                            *local_basis_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),orientation_list(iorientation,:),[1 1 1],[1 1 1]);
                    end
                end
                sol(iorientation,8,3,ix,iy,iz)=summation*gps(ix)*gps(iz);
                %-------------------------------------
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %------------------------------------
                summation=0;
                for ialpha=1:1:3
                    for ibeta=1:1:3
                        term1=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eX,[1 0 0]);
                        term2=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eY,[0 1 0]);
                        term3=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eX,[0 1 0]);
                        term4=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eY,[1 0 0]);
                        determinant=term1*term2-term3*term4;
                        summation=summation+w(ialpha)*w(ibeta)*determinant...
                            *local_basis_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),orientation_list(iorientation,:),[1 1 1],[1 1 1]);
                    end
                end
                sol(iorientation,8,4,ix,iy,iz)=summation*gps(ix)*gps(iy);
                %-------------------------------------
            end
        end
    end
end

%=====================================
%=====================================

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %------------------------------------
                summation=0;
                for ibeta=1:1:3
                    for igamma=1:1:3
                        term1=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 1 0]);
                        term2=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 0 1]);
                        term3=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 0 1]);
                        term4=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 1 0]);
                        determinant=term1*term2-term3*term4;
                        fx=get_fg_dervs_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),orientation_list(iorientation,:),[1 1 1],eX,eY,eZ,1);
                        summation=summation+w(ibeta)*w(igamma)*determinant*fx;
                    end
                end
                sol(iorientation,8,5,ix,iy,iz)=summation*gps(iy)*gps(iz);
                %-------------------------------------
            end
        end
    end
end


for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %------------------------------------
                summation=0;
                for ialpha=1:1:3
                    for igamma=1:1:3
                        term1=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eX,[1 0 0]);
                        term2=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eZ,[0 0 1]);
                        term3=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eX,[0 0 1]);
                        term4=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eZ,[1 0 0]);
                        determinant=term1*term2-term3*term4;
                        fy=get_fg_dervs_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),orientation_list(iorientation,:),[1 1 1],eX,eY,eZ,2);
                        summation=summation+w(ialpha)*w(igamma)*determinant*fy;
                    end
                end
                sol(iorientation,8,6,ix,iy,iz)=summation*gps(ix)*gps(iz);
                %-------------------------------------
            end
        end
    end
end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                %------------------------------------
                summation=0;
                for ialpha=1:1:3
                    for ibeta=1:1:3
                        term1=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eX,[1 0 0]);
                        term2=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eY,[0 1 0]);
                        term3=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eX,[0 1 0]);
                        term4=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),eY,[1 0 0]);
                        determinant=term1*term2-term3*term4;
                        fz=get_fg_dervs_3d(gps_plus(ix,ialpha),gps_plus(iy,ibeta),gps(iz),orientation_list(iorientation,:),[1 1 1],eX,eY,eZ,3);
                        summation=summation+w(ialpha)*w(ibeta)*determinant*fz;
                    end
                end
                sol(iorientation,8,7,ix,iy,iz)=summation*gps(ix)*gps(iy);
                %-------------------------------------
            end
        end
    end
end

end