function sol=type6_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol)

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation1=0;
                summation2=0;
                summation3=0;
                for ialpha=1:1:3
                    for igamma=1:1:3
                        term1=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eX,[1 0 0]);
                        term2=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eZ,[0 0 1]);
                        term3=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eX,[0 0 1]);
                        term4=map_to_global_1_compt_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),eZ,[1 0 0]);
                        determinant=term1*term2-term3*term4;
                        fy=get_fg_dervs_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),...
                            orientation_list(iorientation,:),[1 0 1],eX,eY,eZ,2);
                        summation1=summation1+w(ialpha)*w(igamma)*determinant*fy;
                        fyy=get_fg_dervs_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),...
                            orientation_list(iorientation,:),[1 0 1],eX,eY,eZ,5);
                        summation2=summation2+w(ialpha)*w(igamma)*determinant*fyy;
                        summation3=summation3+w(ialpha)*w(igamma)*determinant...
                            *local_basis_3d(gps_plus(ix,ialpha),gps(iy),gps_plus(iz,igamma),...
                            orientation_list(iorientation,:),[1 0 1],[1 0 1]);
                    end
                end
                mult=gps(ix)*gps(iz);
                sol(iorientation,6,3,ix,iy,iz)=summation1*mult;
                sol(iorientation,6,6,ix,iy,iz)=summation2*mult;
                sol(iorientation,6,1,ix,iy,iz)=summation3*mult;
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
                        orientation_list(iorientation,:),[1 0 1],[1 0 1]);
                    fx=get_fg_dervs_3d(gps(ix),gps(iy),gps_plus(iz,igamma),...
                        orientation_list(iorientation,:),[1 0 1],eX,eY,eZ,1);
                    summation2=summation2+w(igamma)*fx;
                end
                mult=gps(iz)*map_to_global_1_compt_3d(...
                    gps(ix),gps(iy),gps(iz),eZ,[0 0 1]);
                sol(iorientation,6,2,ix,iy,iz)=summation1*mult;
                sol(iorientation,6,5,ix,iy,iz)=summation2*mult;
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
                for ialpha=1:1:3
                    summation1=summation1+w(ialpha)...
                        *local_basis_3d(gps_plus(ix,ialpha),gps(iy),gps(iz),...
                        orientation_list(iorientation,:),[1 0 1],[1 0 1]);
                    fz=get_fg_dervs_3d(gps_plus(ix,ialpha),gps(iy),gps(iz),...
                        orientation_list(iorientation,:),[1 0 1],eX,eY,eZ,3);
                    summation2=summation2+w(ialpha)*fz;
                end
                mult=gps(ix)*map_to_global_1_compt_3d(...
                    gps(ix),gps(iy),gps(iz),eX,[1 0 0]);
                sol(iorientation,6,4,ix,iy,iz)=summation1*mult;
                sol(iorientation,6,7,ix,iy,iz)=summation2*mult;
            end
        end
    end
end

end