function [elements,lbfs1,lbfs2]=analyze_interaction_3d(gbf1,gbf2,nex,ney,nx,ny,nz)

[node1,type1]=analyze_gbs_3d(gbf1);
[node2,type2]=analyze_gbs_3d(gbf2);

[xth1,yth1,zth1]=get_xyzth_3d(node1,nx,ny);
[xth2,yth2,zth2]=get_xyzth_3d(node2,nx,ny);

gbf2_position=positional_element_position(xth1,yth1,zth1,xth2,yth2,zth2,...
    nx,ny,nz);

elements=zeros(1,8);
lbfs1=zeros(1,8);
lbfs2=zeros(1,8);

for iPositional_element=1:1:8
    if gbf2_position(iPositional_element)~=0
        elements(iPositional_element)=get_elements_3d(...
            xth1,yth1,zth1,nex,ney,nx,ny,nz,iPositional_element);
        lbfs1(iPositional_element)=(8-iPositional_element)*8+type1;
        lbfs2(iPositional_element)=(gbf2_position(iPositional_element)...
            -1)*8+type2;
    end
end

% if gbf2_position(1)~=0
%     elements(1)=get_elements_3d(...
%         xth1,yth1,zth1,nex,ney,nx,ny,nz,iPositional_element);
%     lbfs1(1)=(8-iPositional_element)*8+type1;
%     lbfs2(1)=(gbf2_position(iPositional_element)...
%         -1)*8+type2;
% end

end

function positions=positional_element_position(xth1,yth1,zth1,xth2,yth2,zth2,...
    nx,ny,nz)
positions=zeros(8,1);
%=====================================
if zth2<zth1
    if yth2<yth1
        if xth2<xth1
            positions(1)=1;
        elseif xth2==xth1
            if xth2~=1
                positions(1)=2;
            end
            if xth2~=nx
                positions(2)=1;
            end
        elseif xth2>xth1
            positions(2)=2;
        end
    elseif yth2==yth1
        if xth2<xth1
            if yth2~=1
                positions(1)=3;
            end
            if yth2~=ny
                positions(3)=1;
            end
        elseif xth2==xth1
            if yth2~=1
                if xth2~=1
                    positions(1)=4;
                end
                if xth2~=nx
                    positions(2)=3;
                end
            end
            if yth2~=ny
                if xth2~=1
                    positions(3)=2;
                end
                if xth2~=nx
                    positions(4)=1;
                end
            end
        elseif xth2>xth1
            if yth2~=1
                positions(2)=4;
            end
            if yth2~=ny
                positions(4)=2;
            end
        end
	
    elseif yth2>yth1
        if xth2<xth1
            positions(3)=3;
        elseif xth2==xth1
            if xth2~=1
                positions(3)=4;
            end
            if xth2~=nx
                positions(4)=3;
            end
        elseif xth2>xth1
            positions(4)=4;
        end
    end
    %=============================================
elseif zth2==zth1
    if yth2<yth1 %----------------------------------------
        if xth2<xth1
            if zth2~=1
                positions(1)=5;
            end
            if zth2~=nz
                positions(5)=1;
            end
        elseif xth2==xth1
            if zth2~=1
                if xth2~=1
                    positions(1)=6;
                end
                if xth2~=nx
                    positions(2)=5;
                end
            end
            if zth2~=nz
                if xth2~=1
                    positions(5)=2;
                end
                if xth2~=nx
                    positions(6)=1;
                end
            end
        elseif xth2>xth1
            if zth2~=1
                positions(2)=6;
            end
            if zth2~=nz
                positions(6)=2;
            end
        end
    elseif yth2==yth1 %----------------------------------------
        if xth2<xth1
            if zth2~=1
                if yth2~=1
                    positions(1)=7;
                end
                if yth2~=ny
                    positions(3)=5;
                end
            end
            if zth2~=nz
                if yth2~=1
                    positions(5)=3;
                end
                if yth2~=ny
                   positions(7)=1;
                end
            end
        elseif xth2==xth1
            if zth2~=1
                if yth2~=1
                    if xth2~=1
                        positions(1)=8;
                    end
                    if xth2~=nx
                        positions(2)=7;
                    end
                end
                if yth2~=ny
                    if xth2~=1
                        positions(3)=6;
                    end
                    if xth2~=nx
                        positions(4)=5;
                    end
                end
            end
            if zth2~=nz
                if yth2~=1
                    if xth2~=1
                        positions(5)=4;
                    end
                    if xth2~=nx
                        positions(6)=3;
                    end
                end
                if yth2~=ny
                    if xth2~=1
                        positions(7)=2;
                    end
                    if xth2~=nx
                        positions(8)=1;
                    end
                end
            end
        elseif xth2>xth1
            if zth2~=1
                if yth2~=1
                    positions(2)=8;
                end
                if yth2~=ny
                    positions(4)=6;
                end
            end
            if zth2~=nz
                if yth2~=1
                    positions(6)=4;
                end
                if yth2~=ny
                    positions(8)=2;
                end
            end
        end
    elseif yth2>yth1 %----------------------------------------
        if xth2<xth1
            if zth2~=1
                positions(3)=7;
            end
            if zth2~=nz
                positions(7)=3;
            end
        elseif xth2==xth1
            if zth2~=1
                if xth2~=1
                    positions(3)=8;
                end
                if xth2~=nx
                    positions(4)=7;
                end
            end
            if zth2~=nz
                if xth2~=1
                    positions(7)=4;
                end
                if xth2~=nx
                    positions(8)=3;
                end
            end
        elseif xth2>xth1
            if zth2~=1
                positions(4)=8;
            end
            if zth2~=nz
                positions(8)=4;
            end
        end
    end
%================================================================
elseif zth2>zth1
    if yth2<yth1  %----------------------------------------
        if xth2<xth1
            positions(5)=5;
        elseif xth2==xth1
            if xth2~=1
                positions(5)=6;
            end
            if xth2~=nx
                positions(6)=5;
            end
        elseif xth2>xth1
            positions(6)=6;
        end
    elseif yth2==yth1 %----------------------------------------
        if xth2<xth1
            if yth2~=1
                positions(5)=7;
            end
            if yth2~=ny
                positions(7)=5;
            end
        elseif xth2==xth1
            if yth2~=1
                if xth2~=1
                    positions(5)=8;
                end
                if xth2~=nx
                    positions(6)=7;
                end
            end
            if yth2~=ny
                if xth2~=1
                    positions(7)=6;
                end
                if xth2~=nx
                    positions(8)=5;
                end
            end
        elseif xth2>xth1
            if yth2~=1
                positions(6)=8;
            end
            if yth2~=ny
                positions(8)=6;
            end
        end
    elseif yth2>yth1 %----------------------------------------
        if xth2<xth1
            positions(7)=7;
        elseif xth2==xth1
            if xth2~=1
                positions(7)=8;
            end
            if xth2~=nx
                positions(8)=7;
            end
        elseif xth2>xth1
            positions(8)=8;
        end
    end
end
end