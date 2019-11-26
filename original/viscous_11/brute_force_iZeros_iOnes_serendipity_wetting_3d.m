function [iZeros,iOnes,g_bc,bc_wetting]=brute_force_iZeros_iOnes_serendipity_wetting_3d(...
    gbfs1,gbfs2,nnz_,nworkers,nx,ny,nz,...
    g_x_minus,g_y_minus,g_z_minus,g_x_plus,g_y_plus,g_z_plus,...
    h_x_minus,h_y_minus,h_z_minus,h_x_plus,h_y_plus,h_z_plus)

%This one has very ugly algorithm


allocation=wise_task_splitter(nnz_,nworkers);

[~,n]=size(allocation);
extra=mod(nnz_,nworkers);
load_less=(nnz_-extra)/nworkers;
if extra==0
    load_more=load_less;
else
    load_more=load_less+1;
end
temp=zeros(nworkers,n);

parfor worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    else
        load=load_less;
    end
    temp(worker,:)=get_keyVals_brute_force_wetting(...
        n,allocation(worker,:),load,gbfs1,gbfs2,nx,ny,nz);
end
    


iZeros=zeros(1,nnz_);
iOnes=zeros(1,4*nx*ny*nz);
i_iZeros=0;
i_iOnes=0;
index=0;

i_g_xdir_yth_1=0;
i_g_xdir_zth_1=1;

i_g_xdir_yth_2=0;
i_g_xdir_zth_2=1;

i_g_ydir_xth_1=0;
i_g_ydir_zth_1=1;

i_g_ydir_xth_2=0;
i_g_ydir_zth_2=1;

i_g_zdir_xth_1=0;
i_g_zdir_yth_1=1;

i_g_zdir_xth_2=0;
i_g_zdir_yth_2=1;


g_bc=zeros(2,nnz_);
i_g_bc=0;

bc_wetting=zeros(4,nnz_);

for worker=1:1:nworkers
    if worker<=extra
        load=load_more;
    else
        load=load_less;
    end
    for i=1:1:load
        index=index+1;
        if temp(worker,i)==1
            i_iZeros=i_iZeros+1;
            iZeros(i_iZeros)=index;
        elseif temp(worker,i)==2
            i_iOnes=i_iOnes+1;
            iOnes(i_iOnes)=index;
        elseif temp(worker,i)==3
            i_g_bc=i_g_bc+1;
            i_g_xdir_yth_1=i_g_xdir_yth_1+1;
            g_bc(1,i_g_bc)=index;
            g_bc(2,i_g_bc)=-g_x_minus(i_g_xdir_yth_1,i_g_xdir_zth_1);
            
            bc_wetting(1,i_g_bc)=gbfs1(index);
            bc_wetting(2,i_g_bc)=gbfs1(index)-1;
            bc_wetting(3,i_g_bc)=h_x_minus(i_g_xdir_yth_1,i_g_xdir_zth_1);
            bc_wetting(4,i_g_bc)=g_x_minus(i_g_xdir_yth_1,i_g_xdir_zth_1);
            
            if i_g_xdir_yth_1==ny
                i_g_xdir_yth_1=0;
                i_g_xdir_zth_1=i_g_xdir_zth_1+1;
                
            end
            
        elseif temp(worker,i)==4
            i_g_bc=i_g_bc+1;
            i_g_xdir_yth_2=i_g_xdir_yth_2+1;
            g_bc(1,i_g_bc)=index;
            g_bc(2,i_g_bc)=-g_x_plus(i_g_xdir_yth_2,i_g_xdir_zth_2);
            
            bc_wetting(1,i_g_bc)=gbfs1(index);
            bc_wetting(2,i_g_bc)=gbfs1(index)-1;
            bc_wetting(3,i_g_bc)=h_x_plus(i_g_xdir_yth_2,i_g_xdir_zth_2);
            bc_wetting(4,i_g_bc)=g_x_plus(i_g_xdir_yth_2,i_g_xdir_zth_2);
            
            if i_g_xdir_yth_2==ny
                i_g_xdir_yth_2=0;
                i_g_xdir_zth_2=i_g_xdir_zth_2+1;
                
            end
            
        elseif temp(worker,i)==5
            i_g_bc=i_g_bc+1;
            i_g_ydir_xth_1=i_g_ydir_xth_1+1;
            g_bc(1,i_g_bc)=index;
            g_bc(2,i_g_bc)=-g_y_minus(i_g_ydir_xth_1,i_g_ydir_zth_1);
            
            bc_wetting(1,i_g_bc)=gbfs1(index);
            bc_wetting(2,i_g_bc)=gbfs1(index)-2;
            bc_wetting(3,i_g_bc)=h_y_minus(i_g_ydir_xth_1,i_g_ydir_zth_1);
            bc_wetting(4,i_g_bc)=g_y_minus(i_g_ydir_xth_1,i_g_ydir_zth_1);
            
            if i_g_ydir_xth_1==nx
                i_g_ydir_xth_1=0;
                i_g_ydir_zth_1=i_g_ydir_zth_1+1;
                
            end
            
        elseif temp(worker,i)==6
            i_g_bc=i_g_bc+1;
            i_g_ydir_xth_2=i_g_ydir_xth_2+1;
            g_bc(1,i_g_bc)=index;
            g_bc(2,i_g_bc)=-g_y_plus(i_g_ydir_xth_2,i_g_ydir_zth_2);
            
            bc_wetting(1,i_g_bc)=gbfs1(index);
            bc_wetting(2,i_g_bc)=gbfs1(index)-2;
            bc_wetting(3,i_g_bc)=h_y_plus(i_g_ydir_xth_2,i_g_ydir_zth_2);
            bc_wetting(4,i_g_bc)=g_y_plus(i_g_ydir_xth_2,i_g_ydir_zth_2);
            
            if i_g_ydir_xth_2==nx
                i_g_ydir_xth_2=0;
                i_g_ydir_zth_2=i_g_ydir_zth_2+1;
                
            end
            
        elseif temp(worker,i)==7
            i_g_bc=i_g_bc+1;
            i_g_zdir_xth_1=i_g_zdir_xth_1+1;
            g_bc(1,i_g_bc)=index;
            g_bc(2,i_g_bc)=-g_z_minus(i_g_zdir_xth_1,i_g_zdir_yth_1);
            
            bc_wetting(1,i_g_bc)=gbfs1(index);
            bc_wetting(2,i_g_bc)=gbfs1(index)-3;
            bc_wetting(3,i_g_bc)=h_z_minus(i_g_zdir_xth_1,i_g_zdir_yth_1);
            bc_wetting(4,i_g_bc)=g_z_minus(i_g_zdir_xth_1,i_g_zdir_yth_1);
            
            if i_g_zdir_xth_1==nx
                i_g_zdir_xth_1=0;
                i_g_zdir_yth_1=i_g_zdir_yth_1+1;
                
            end
            
        elseif temp(worker,i)==8
            i_g_bc=i_g_bc+1;
            i_g_zdir_xth_2=i_g_zdir_xth_2+1;
            g_bc(1,i_g_bc)=index;
            g_bc(2,i_g_bc)=-g_z_plus(i_g_zdir_xth_2,i_g_zdir_yth_2);
            
            bc_wetting(1,i_g_bc)=gbfs1(index);
            bc_wetting(2,i_g_bc)=gbfs1(index)-3;
            bc_wetting(3,i_g_bc)=h_z_plus(i_g_zdir_xth_2,i_g_zdir_yth_2);
            bc_wetting(4,i_g_bc)=g_z_plus(i_g_zdir_xth_2,i_g_zdir_yth_2);
            
            if i_g_zdir_xth_2==nx
                i_g_zdir_xth_2=0;
                i_g_zdir_yth_2=i_g_zdir_yth_2+1;
                
            end
            
        end
    end
end
reducer=nnz(iZeros)+1;
iZeros(reducer:1:nnz_)=[];
reducer=nnz(iOnes)+1;
iOnes(reducer:1:4*nx*ny*nz)=[];
reducer=nnz(g_bc(1,:));
g_bc=g_bc(:,1:1:reducer);
bc_wetting=bc_wetting(:,1:1:reducer);

end

function sol=get_keyVals_brute_force_wetting(n,keys,load,gbfs1,gbfs2,nx,ny,nz)

sol=zeros(1,n);
for i=1:1:load
    [node1,type1]=analyze_gbs_serendipity(gbfs1(keys(i)));
    [xth1,yth1,zth1]=get_xyzth_3d(node1,nx,ny);
    if type1==2 && (xth1==1 || xth1==nx)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            [node2,type2]=analyze_gbs_serendipity(gbfs2(keys(i)));
            if node2==node1 && type2==1
                if xth1==1
                    sol(i)=3;
                elseif xth1==nx
                    sol(i)=4;
                end
            else
                sol(i)=1;
            end
        end

    
    elseif type1==3 && (yth1==1 || yth1==ny)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            [node2,type2]=analyze_gbs_serendipity(gbfs2(keys(i)));
            if node2==node1 && type2==1
                if yth1==1
                    sol(i)=5;
                elseif yth1==ny
                    sol(i)=6;
                end
            else
                sol(i)=1;
            end
        end
            

    elseif type1==4 && (zth1==1 || zth1==nz)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            [node2,type2]=analyze_gbs_serendipity(gbfs2(keys(i)));
            if node2==node1 && type2==1
                if zth1==1
                    sol(i)=7;
                elseif zth1==nz
                    sol(i)=8;
                end
            else
                sol(i)=1;
            end
        end
    end
    
end
end

function sol=get_keyVals_brute_force(n,keys,load,gbfs1,gbfs2,nx,ny,nz)

sol=zeros(1,n);
for i=1:1:load
    [node1,type1]=analyze_gbs_serendipity(gbfs1(keys(i)));
    [xth1,yth1,zth1]=get_xyzth_3d(node1,nx,ny);
    if type1==2 && (xth1==1 || xth1==nx)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            sol(i)=1;
        end

    
    elseif type1==3 && (yth1==1 || yth1==ny)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            sol(i)=1;
        end
            

    elseif type1==4 && (zth1==1 || zth1==nz)

        if gbfs1(keys(i))==gbfs2(keys(i))
            sol(i)=2;
        else
            sol(i)=1;
        end

    end
end
end