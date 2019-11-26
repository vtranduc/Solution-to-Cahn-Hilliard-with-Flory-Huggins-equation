function sol=type1_transformer_3d(gps,eX,eY,eZ,orientation_list,sol)

for iorientation=1:1:8
    for iz=1:1:3
        for iy=1:1:3
            for ix=1:1:3
                sol(iorientation,1,1,ix,iy,iz)=...
                    local_basis_3d(gps(ix),gps(iy),gps(iz),...
                    orientation_list(iorientation,:),[0 0 0],[0 0 0]);
                df=get_fg_dervs_3d(gps(ix),gps(iy),gps(iz),...
                    orientation_list(iorientation,:),[0 0 0],eX,eY,eZ,1:1:6);
                for idf=1:1:6
                    sol(iorientation,1,idf+1,ix,iy,iz)=df(idf);
                end
            end
        end
    end
end

end