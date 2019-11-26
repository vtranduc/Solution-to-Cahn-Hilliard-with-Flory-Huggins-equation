function [elements,lbfs1,lbfs2,mx1,mx2,mx3,mx4,mx5,mx6]=...
    analyze_interaction_sph(gbf1,gbf2,ne,nex,ney,nez,n,nx,ny,nz,neight,...
    nexney,n_eSurf,n_nSurf,nx_)

if gbf1<=neight && gbf2<=neight
    [elements,lbfs1,lbfs2,xth1,yth1,zth1,xth2,yth2,zth2,...
        gbf2_position,type1,type2]=...
        analyze_interaction_traditional(gbf1,gbf2,nex,ney,nx,ny,nz);
    
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    %------------------------------------------------------
    
%     if (xth1==1 && xth2==1) || (xth1==nx && xth2==nx) || (yth1==1 && yth2==1) || (yth1==ny && yth2==ny) || (zth1==1 && zth2==1) || (zth1==nz && zth2==nz)
%         pos_mx=zeros(8,2);
%         for i=1:1:8
%             if gbf2_position(i)~=0
%                 pos_mx(i,1)=9-i;
%                 pos_mx(i,2)=gbf2_position(i);
%             end
%         end
%     else
%         mx1=NaN;mx2=NaN;mx3=NaN;mx4=NaN;mx5=NaN;mx6=NaN;
%         return
%     end
    
    %--------------------
    if xth1==1 && xth2==1
        mx1=zeros(4,3);
        index=0;
        for iPos=[2 4 6 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx1(index,1)=esurfs(1)+ne;
                mx1(index,2)=mirror(9-iPos,1);
                mx1(index,3)=mirror(gbf2_position(iPos),1);
                mx1(index,2)=8*(mx1(index,2)-1)+type1;
                mx1(index,3)=8*(mx1(index,3)-1)+type2;
            end
        end
    else
        mx1=NaN;
    end
    %--------------------
    
    %--------------------
    if xth1==nx && xth2==nx
        mx2=zeros(4,3);
        index=0;
        for iPos=[1 3 5 7]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx2(index,1)=esurfs(1)+ne;
                mx2(index,2)=mirror(9-iPos,2);
                mx2(index,3)=mirror(gbf2_position(iPos),2);
                mx2(index,2)=8*(mx2(index,2)-1)+type1;
                mx2(index,3)=8*(mx2(index,3)-1)+type2;
            end
        end
    else
        mx2=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==1 && yth2==1
        mx3=zeros(4,3);
        index=0;
        for iPos=[3 4 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx3(index,1)=esurfs(2)+ne;
                mx3(index,2)=mirror(9-iPos,3);
                mx3(index,3)=mirror(gbf2_position(iPos),3);
                mx3(index,2)=8*(mx3(index,2)-1)+type1;
                mx3(index,3)=8*(mx3(index,3)-1)+type2;
            end
        end
    else
        mx3=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==ny && yth2==ny
        mx4=zeros(4,3);
        index=0;
        for iPos=[1 2 5 6]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx4(index,1)=esurfs(2)+ne;
                mx4(index,2)=mirror(9-iPos,4);
                mx4(index,3)=mirror(gbf2_position(iPos),4);
                mx4(index,2)=8*(mx4(index,2)-1)+type1;
                mx4(index,3)=8*(mx4(index,3)-1)+type2;
            end
        end
    else
        mx4=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==1 && zth2==1
        mx5=zeros(4,3);
        index=0;
        for iPos=[5 6 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx5(index,1)=esurfs(3)+ne;
                mx5(index,2)=mirror(9-iPos,5);
                mx5(index,3)=mirror(gbf2_position(iPos),5);
                mx5(index,2)=8*(mx5(index,2)-1)+type1;
                mx5(index,3)=8*(mx5(index,3)-1)+type2;
            end
        end
    else
        mx5=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==nz && zth2==nz
        mx6=zeros(4,3);
        index=0;
        for iPos=[1 2 3 4]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx6(index,1)=esurfs(3)+ne;
                mx6(index,2)=mirror(9-iPos,6);
                mx6(index,3)=mirror(gbf2_position(iPos),6);
                mx6(index,2)=8*(mx6(index,2)-1)+type1;
                mx6(index,3)=8*(mx6(index,3)-1)+type2;
            end
        end
    else
        mx6=NaN;
    end
    %--------------------
    
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    
elseif gbf1>neight && gbf2<=neight
    
    [node1,type1]=analyze_gbs_3d(gbf1);
    [xth1,yth1,zth1]=get_root_of_extrusion_sph(node1,n_nSurf,n,nx,ny,nz,nx_);
    gbf0=(get_node_3d(xth1,yth1,zth1,nx,ny)-1)*8+type1;
    [elements,~,~,xth1,yth1,zth1,xth2,yth2,zth2,...
        gbf2_position,type1,type2]=...
        analyze_interaction_traditional(gbf0,gbf2,nex,ney,nx,ny,nz);
    
    %--------------------
    if xth1==1 && xth2==1
        mx1=zeros(4,3);
        index=0;
        for iPos=[2 4 6 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx1(index,1)=esurfs(1)+ne;
                mx1(index,2)=9-iPos;
                mx1(index,3)=mirror(gbf2_position(iPos),1);
                mx1(index,2)=8*(mx1(index,2)-1)+type1;
                mx1(index,3)=8*(mx1(index,3)-1)+type2;
            end
        end
    else
        mx1=NaN;
    end
    %--------------------
    
    %--------------------
    if xth1==nx && xth2==nx
        mx2=zeros(4,3);
        index=0;
        for iPos=[1 3 5 7]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx2(index,1)=esurfs(1)+ne;
                mx2(index,2)=9-iPos;
                mx2(index,3)=mirror(gbf2_position(iPos),2);
                mx2(index,2)=8*(mx2(index,2)-1)+type1;
                mx2(index,3)=8*(mx2(index,3)-1)+type2;
            end
        end
    else
        mx2=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==1 && yth2==1
        mx3=zeros(4,3);
        index=0;
        for iPos=[3 4 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx3(index,1)=esurfs(2)+ne;
                mx3(index,2)=9-iPos;
                mx3(index,3)=mirror(gbf2_position(iPos),3);
                mx3(index,2)=8*(mx3(index,2)-1)+type1;
                mx3(index,3)=8*(mx3(index,3)-1)+type2;
            end
        end
    else
        mx3=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==ny && yth2==ny
        mx4=zeros(4,3);
        index=0;
        for iPos=[1 2 5 6]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx4(index,1)=esurfs(2)+ne;
                mx4(index,2)=9-iPos;
                mx4(index,3)=mirror(gbf2_position(iPos),4);
                mx4(index,2)=8*(mx4(index,2)-1)+type1;
                mx4(index,3)=8*(mx4(index,3)-1)+type2;
            end
        end
    else
        mx4=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==1 && zth2==1
        mx5=zeros(4,3);
        index=0;
        for iPos=[5 6 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx5(index,1)=esurfs(3)+ne;
                mx5(index,2)=9-iPos;
                mx5(index,3)=mirror(gbf2_position(iPos),5);
                mx5(index,2)=8*(mx5(index,2)-1)+type1;
                mx5(index,3)=8*(mx5(index,3)-1)+type2;
            end
        end
    else
        mx5=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==nz && zth2==nz
        mx6=zeros(4,3);
        index=0;
        for iPos=[1 2 3 4]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx6(index,1)=esurfs(3)+ne;
                mx6(index,2)=9-iPos;
                mx6(index,3)=mirror(gbf2_position(iPos),6);
                mx6(index,2)=8*(mx6(index,2)-1)+type1;
                mx6(index,3)=8*(mx6(index,3)-1)+type2;
            end
        end
    else
        mx6=NaN;
    end
    %--------------------
    
    elements=NaN;lbfs1=NaN;lbfs2=NaN;
    
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    
elseif gbf1<=neight && gbf2>neight
    
    [node2,type2]=analyze_gbs_3d(gbf2);
    [xth2,yth2,zth2]=get_root_of_extrusion_sph(node2,n_nSurf,n,nx,ny,nz,nx_);
    gbf0=(get_node_3d(xth2,yth2,zth2,nx,ny)-1)*8+type2;
    [elements,~,~,xth1,yth1,zth1,xth2,yth2,zth2,...
        gbf2_position,type1,type2]=...
        analyze_interaction_traditional(gbf1,gbf0,nex,ney,nx,ny,nz);
    
    %--------------------
    if xth1==1 && xth2==1
        mx1=zeros(4,3);
        index=0;
        for iPos=[2 4 6 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx1(index,1)=esurfs(1)+ne;
                mx1(index,2)=mirror(9-iPos,1);
                mx1(index,3)=gbf2_position(iPos);
                mx1(index,2)=8*(mx1(index,2)-1)+type1;
                mx1(index,3)=8*(mx1(index,3)-1)+type2;
            end
        end
    else
        mx1=NaN;
    end
    %--------------------
    
    %--------------------
    if xth1==nx && xth2==nx
        mx2=zeros(4,3);
        index=0;
        for iPos=[1 3 5 7]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx2(index,1)=esurfs(1)+ne;
                mx2(index,2)=mirror(9-iPos,2);
                mx2(index,3)=gbf2_position(iPos);
                mx2(index,2)=8*(mx2(index,2)-1)+type1;
                mx2(index,3)=8*(mx2(index,3)-1)+type2;
            end
        end
    else
        mx2=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==1 && yth2==1
        mx3=zeros(4,3);
        index=0;
        for iPos=[3 4 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx3(index,1)=esurfs(2)+ne;
                mx3(index,2)=mirror(9-iPos,3);
                mx3(index,3)=gbf2_position(iPos);
                mx3(index,2)=8*(mx3(index,2)-1)+type1;
                mx3(index,3)=8*(mx3(index,3)-1)+type2;
            end
        end
    else
        mx3=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==ny && yth2==ny
        mx4=zeros(4,3);
        index=0;
        for iPos=[1 2 5 6]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx4(index,1)=esurfs(2)+ne;
                mx4(index,2)=mirror(9-iPos,4);
                mx4(index,3)=gbf2_position(iPos);
                mx4(index,2)=8*(mx4(index,2)-1)+type1;
                mx4(index,3)=8*(mx4(index,3)-1)+type2;
            end
        end
    else
        mx4=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==1 && zth2==1
        mx5=zeros(4,3);
        index=0;
        for iPos=[5 6 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx5(index,1)=esurfs(3)+ne;
                mx5(index,2)=mirror(9-iPos,5);
                mx5(index,3)=gbf2_position(iPos);
                mx5(index,2)=8*(mx5(index,2)-1)+type1;
                mx5(index,3)=8*(mx5(index,3)-1)+type2;
            end
        end
    else
        mx5=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==nz && zth2==nz
        mx6=zeros(4,3);
        index=0;
        for iPos=[1 2 3 4]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx6(index,1)=esurfs(3)+ne;
                mx6(index,2)=mirror(9-iPos,6);
                mx6(index,3)=gbf2_position(iPos);
                mx6(index,2)=8*(mx6(index,2)-1)+type1;
                mx6(index,3)=8*(mx6(index,3)-1)+type2;
            end
        end
    else
        mx6=NaN;
    end
    %--------------------
    
    elements=NaN;lbfs1=NaN;lbfs2=NaN;
    
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    %===================================================================
    
elseif gbf1>neight && gbf2>neight
    
    [node1,type1]=analyze_gbs_3d(gbf1);
    [xth1,yth1,zth1]=get_root_of_extrusion_sph(node1,n_nSurf,n,nx,ny,nz,nx_);
    gbf1_=(get_node_3d(xth1,yth1,zth1,nx,ny)-1)*8+type1;
    
    [node2,type2]=analyze_gbs_3d(gbf2);
    [xth2,yth2,zth2]=get_root_of_extrusion_sph(node2,n_nSurf,n,nx,ny,nz,nx_);
    gbf2_=(get_node_3d(xth2,yth2,zth2,nx,ny)-1)*8+type2;
    
    [elements,~,~,xth1,yth1,zth1,xth2,yth2,zth2,...
        gbf2_position,type1,type2]=...
        analyze_interaction_traditional(gbf1_,gbf2_,nex,ney,nx,ny,nz);
    
    %--------------------
    if xth1==1 && xth2==1
        mx1=zeros(4,3);
        index=0;
        for iPos=[2 4 6 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx1(index,1)=esurfs(1)+ne;
                mx1(index,2)=9-iPos;
                mx1(index,3)=gbf2_position(iPos);
                mx1(index,2)=8*(mx1(index,2)-1)+type1;
                mx1(index,3)=8*(mx1(index,3)-1)+type2;
            end
        end
    else
        mx1=NaN;
    end
    %--------------------
    
    %--------------------
    if xth1==nx && xth2==nx
        mx2=zeros(4,3);
        index=0;
        for iPos=[1 3 5 7]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx2(index,1)=esurfs(1)+ne;
                mx2(index,2)=9-iPos;
                mx2(index,3)=gbf2_position(iPos);
                mx2(index,2)=8*(mx2(index,2)-1)+type1;
                mx2(index,3)=8*(mx2(index,3)-1)+type2;
            end
        end
    else
        mx2=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==1 && yth2==1
        mx3=zeros(4,3);
        index=0;
        for iPos=[3 4 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx3(index,1)=esurfs(2)+ne;
                mx3(index,2)=9-iPos;
                mx3(index,3)=gbf2_position(iPos);
                mx3(index,2)=8*(mx3(index,2)-1)+type1;
                mx3(index,3)=8*(mx3(index,3)-1)+type2;
            end
        end
    else
        mx3=NaN;
    end
    %--------------------
    
    %--------------------
    if yth1==ny && yth2==ny
        mx4=zeros(4,3);
        index=0;
        for iPos=[1 2 5 6]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx4(index,1)=esurfs(2)+ne;
                mx4(index,2)=9-iPos;
                mx4(index,3)=gbf2_position(iPos);
                mx4(index,2)=8*(mx4(index,2)-1)+type1;
                mx4(index,3)=8*(mx4(index,3)-1)+type2;
            end
        end
    else
        mx4=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==1 && zth2==1
        mx5=zeros(4,3);
        index=0;
        for iPos=[5 6 7 8]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx5(index,1)=esurfs(3)+ne;
                mx5(index,2)=9-iPos;
                mx5(index,3)=gbf2_position(iPos);
                mx5(index,2)=8*(mx5(index,2)-1)+type1;
                mx5(index,3)=8*(mx5(index,3)-1)+type2;
            end
        end
    else
        mx5=NaN;
    end
    %--------------------
    
    %--------------------
    if zth1==nz && zth2==nz
        mx6=zeros(4,3);
        index=0;
        for iPos=[1 2 3 4]
            index=index+1;
            if elements(iPos)~=0
                esurfs=get_extruding_element_sph(elements(iPos),nex,ney,nez,nexney,n_eSurf);
                mx6(index,1)=esurfs(3)+ne;
                mx6(index,2)=9-iPos;
                mx6(index,3)=gbf2_position(iPos);
                mx6(index,2)=8*(mx6(index,2)-1)+type1;
                mx6(index,3)=8*(mx6(index,3)-1)+type2;
            end
        end
    else
        mx6=NaN;
    end
    %--------------------
    
    elements=NaN;lbfs1=NaN;lbfs2=NaN;
    
end


end


function [elements,lbfs1,lbfs2,xth1,yth1,zth1,xth2,yth2,zth2,...
    gbf2_position,type1,type2]=...
    analyze_interaction_traditional(gbf1,gbf2,nex,ney,nx,ny,nz)

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