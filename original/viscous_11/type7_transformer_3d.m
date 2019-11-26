function sol=type7_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol)

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation1=0;
                summation2=0;
                summation3=0;
                for ibeta=1:1:3
                    for igamma=1:1:3
                        term1=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 1 0]);
                        term2=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 0 1]);
                        term3=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eY,[0 0 1]);
                        term4=map_to_global_1_compt_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),eZ,[0 1 0]);
                        determinant=term1*term2-term3*term4;
                        fx=get_fg_dervs_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),...
                            orientation_list(iorientation,:),[0 1 1],eX,eY,eZ,1);
                        summation1=summation1+w(ibeta)*w(igamma)*determinant*fx;
                        fxx=get_fg_dervs_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),...
                            orientation_list(iorientation,:),[0 1 1],eX,eY,eZ,4);
                        summation2=summation2+w(ibeta)*w(igamma)*determinant*fxx;
                        summation3=summation3+w(ibeta)*w(igamma)*determinant...
                            *local_basis_3d(gps(ix),gps_plus(iy,ibeta),gps_plus(iz,igamma),...
                            orientation_list(iorientation,:),[0 1 1],[0 1 1]);
                    end
                end
                mult=gps(iy)*gps(iz);
                sol(iorientation,7,2,ix,iy,iz)=summation1*mult;
                sol(iorientation,7,5,ix,iy,iz)=summation2*mult;
                sol(iorientation,7,1,ix,iy,iz)=summation3*mult;
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
                for igamma=1:1:3
                    summation1=summation1+w(igamma)...
                        *local_basis_3d(gps(ix),gps(iy),gps_plus(iz,igamma),...
                        orientation_list(iorientation,:),[0 1 1],[0 1 1]);
                    fy=get_fg_dervs_3d(gps(ix),gps(iy),gps_plus(iz,igamma),...
                        orientation_list(iorientation,:),[0 1 1],eX,eY,eZ,2);
                    summation2=summation2+w(igamma)*fy;
                end
                mult=gps(iz)*map_to_global_1_compt_3d(...
                    gps(ix),gps(iy),gps(iz),eZ,[0 0 1]);
                sol(iorientation,7,3,ix,iy,iz)=summation1*mult;
                sol(iorientation,7,6,ix,iy,iz)=summation2*mult;
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
                for ibeta=1:1:3
                    summation1=summation1+w(ibeta)...
                        *local_basis_3d(gps(ix),gps_plus(iy,ibeta),gps(iz),...
                        orientation_list(iorientation,:),[0 1 1],[0 1 1]);
                    fz=get_fg_dervs_3d(gps(ix),gps_plus(iy,ibeta),gps(iz),...
                        orientation_list(iorientation,:),[0 1 1],eX,eY,eZ,3);
                    summation2=summation2+w(ibeta)*fz;
                end
                mult=gps(iy)*map_to_global_1_compt_3d(...
                    gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
                sol(iorientation,7,4,ix,iy,iz)=summation1*mult;
                sol(iorientation,7,7,ix,iy,iz)=summation2*mult;
            end
        end
    end
end

end