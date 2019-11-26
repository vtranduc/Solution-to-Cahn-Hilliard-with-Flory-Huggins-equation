function sol=type3_transformer_3d(gps,gps_plus,w,eX,eY,eZ,orientation_list,sol)

 for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                summation1=0;
                summation2=0;
                summation3=0;
                summation4=0;
                summation5=0;
                for ibeta=1:1:3
                    summation1=summation1+w(ibeta)...
                        *local_basis_3d(gps(ix),gps_plus(iy,ibeta),gps(iz),...
                        orientation_list(iorientation,:),[0 1 0],[0 1 0]);
                    df=get_fg_dervs_3d(gps(ix),gps_plus(iy,ibeta),gps(iz),...
                        orientation_list(iorientation,:),[0 1 0],eX,eY,eZ,[1 3 4 6]);
                    summation2=summation2+w(ibeta)*df(1);
                    summation3=summation3+w(ibeta)*df(2);
                    summation4=summation4+w(ibeta)*df(3);
                    summation5=summation5+w(ibeta)*df(4);
                end
                mult=gps(iy)*map_to_global_1_compt_3d(...
                    gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
                sol(iorientation,3,1,ix,iy,iz)=summation1*mult;
                sol(iorientation,3,2,ix,iy,iz)=summation2*mult;
                sol(iorientation,3,4,ix,iy,iz)=summation3*mult;
                sol(iorientation,3,5,ix,iy,iz)=summation4*mult;
                sol(iorientation,3,7,ix,iy,iz)=summation5*mult;
            end
        end
    end
 end

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                sol(iorientation,3,3,ix,iy,iz)=local_basis_3d(gps(ix),gps(iy),gps(iz),...
                    orientation_list(iorientation,:),[0 1 0],[0 1 0]);
                sol(iorientation,3,6,ix,iy,iz)=get_fg_dervs_3d(gps(ix),gps(iy),gps(iz),...
                    orientation_list(iorientation,:),[0 1 0],eX,eY,eZ,2);
            end
        end
    end
end

end