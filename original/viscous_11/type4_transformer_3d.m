function sol=type4_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol)

 for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation1=0;
                summation2=0;
                summation3=0;
                summation4=0;
                summation5=0;
                for igamma=1:1:3
                    summation1=summation1+w(igamma)...
                        *local_basis_3d(gps(ix),gps(iy),gps_plus(iz,igamma),...
                        orientation_list(iorientation,:),[0 0 1],[0 0 1]);
                    df=get_fg_dervs_3d(gps(ix),gps(iy),gps_plus(iz,igamma),...
                        orientation_list(iorientation,:),[0 0 1],eX,eY,eZ,[1 2 4 5]);
                    summation2=summation2+w(igamma)*df(1);
                    summation3=summation3+w(igamma)*df(2);
                    summation4=summation4+w(igamma)*df(3);
                    summation5=summation5+w(igamma)*df(4);
                end
                mult=gps(iz)*map_to_global_1_compt_3d(...
                    gps(ix),gps(iy),gps(iz),eZ,[0 0 1]);
                sol(iorientation,4,1,ix,iy,iz)=summation1*mult;
                sol(iorientation,4,2,ix,iy,iz)=summation2*mult;
                sol(iorientation,4,3,ix,iy,iz)=summation3*mult;
                sol(iorientation,4,5,ix,iy,iz)=summation4*mult;
                sol(iorientation,4,6,ix,iy,iz)=summation5*mult;
            end
        end
    end
 end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                
                sol(iorientation,4,4,ix,iy,iz)=local_basis_3d(gps(ix),gps(iy),gps(iz),...
                    orientation_list(iorientation,:),[0 0 1],[0 0 1]);
                
                
                sol(iorientation,4,7,ix,iy,iz)=get_fg_dervs_3d(gps(ix),gps(iy),gps(iz),...
                    orientation_list(iorientation,:),[0 0 1],eX,eY,eZ,3);
            end
        end
    end
end

end