function [diffT,chi]=create_temperature_gradient_3d(parallel_computing,shape,T,diff,entropy,T_theta,nx,ny,nexney,ne,nex,X,Y,Z)
%This will only work assuming the meshes have been well defined

if length(T)==1
    diffT=diff*T;
    chi=0.5-entropy*(1-T_theta/T);
else
    gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
    diffT=zeros(ne,3,3,3);
    chi=zeros(ne,3,3,3);
    if parallel_computing==1
        if shape==1
            slope=T(2)-T(1);
            x_dist=X(2)-X(1);
            gps=gps*x_dist;
            start=T(1);
            parfor e=1:1:ne
                [temp_pts,chi_pts]=horizontal_gradient_3d(e,nex,gps,x_dist,start,slope,entropy,T_theta);
                diffT(e,:,:,:)=temp_pts*diff;
                chi(e,:,:,:)=chi_pts;
            end
        elseif shape==2 || shape==3
            start=T(1);
            slope=T(2)-T(1);
            parfor e=1:1:ne
                [temp_pts,chi_pts]=radial_gradient_3d(e,nex,nx,ny,nexney,X,Y,Z,start,slope,gps,entropy,T_theta);
                diffT(e,:,:,:)=temp_pts*diff;
                chi(e,:,:,:)=chi_pts;
            end
        else
            error('Only 2 shapes, 1 (cube) and 2 (sphere) are available')
        end
    elseif parallel_computing==0
        if shape==1
            slope=T(2)-T(1);
            x_dist=X(2)-X(1);
            gps=gps*x_dist;
            for e=1:1:ne
                n1xth=get_n1xth_3d(e,nex);
                n1dist=(n1xth-1)*x_dist;
                for ix=1:1:3
                    for iy=1:1:3
                        for iz=1:1:3
                            diffT(e,ix,iy,iz)=T(1)+slope*(n1dist+gps(ix));
                            chi(e,ix,iy,iz)=0.5-entropy*(1-T_theta/diffT(e,ix,iy,iz));
                        end
                    end
                end
            end
            diffT=diffT*diff;
        elseif shape==2 || shape==3
            start=T(1);
            slope=T(2)-T(1);
            for e=1:1:ne
                [temp_pts,chi_pts]=radial_gradient_3d(e,nex,nx,ny,nexney,X,Y,Z,start,slope,gps,entropy,T_theta);
                diffT(e,:,:,:)=temp_pts*diff;
                chi(e,:,:,:)=chi_pts;
            end
        else
            error('Only 2 shapes, 1 (cube) and 2 (sphere) are available')
        end
    end
end

end

function [temp_pts,chi_pts]=horizontal_gradient_3d(e,nex,gps,x_dist,start,slope,entropy,T_theta)
temp_pts=zeros(3,3,3);
chi_pts=zeros(3,3,3);
n1xth=get_n1xth_3d(e,nex);
n1dist=(n1xth-1)*x_dist;
for ix=1:1:3
    for iy=1:1:3
        for iz=1:1:3
            temp_pts(ix,iy,iz)=start+slope*(n1dist+gps(ix));
            chi_pts(ix,iy,iz)=0.5-entropy*(1-T_theta/temp_pts(ix,iy,iz));
        end
    end
end
end

function [temp_pts,chi_pts]=radial_gradient_3d(e,nex,nx,ny,nexney,X,Y,Z,start,slope,gps,entropy,T_theta)
temp_pts=zeros(3,3,3);
chi_pts=zeros(3,3,3);
[eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z);
for igamma=1:1:3
    for ibeta=1:1:3
        for ialpha=1:1:3
            x=map_to_global_1_compt_3d(gps(ialpha),gps(ibeta),gps(igamma),eX,[0 0 0]);
            y=map_to_global_1_compt_3d(gps(ialpha),gps(ibeta),gps(igamma),eY,[0 0 0]);
            z=map_to_global_1_compt_3d(gps(ialpha),gps(ibeta),gps(igamma),eZ,[0 0 0]);
            radius=sqrt(x^2+y^2+z^2);
            temp_pts(ialpha,ibeta,igamma)=start+slope*radius;
            chi_pts(ialpha,ibeta,igamma)=0.5-entropy*(1-T_theta/temp_pts(ialpha,ibeta,igamma));
        end
    end
end

end