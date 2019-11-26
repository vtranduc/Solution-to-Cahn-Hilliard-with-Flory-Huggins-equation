function sol=compute_weights_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list)

sol=zeros(8,8,7);
for iorientation=1:1:8
    transformation=type1_transformer_specific_3d(alpha,beta,gamma,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,1,iorder)=transformation(iorder);
    end
    transformation=type2_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,2,iorder)=transformation(iorder);
    end
    transformation=type3_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,3,iorder)=transformation(iorder);
    end
    transformation=type4_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,4,iorder)=transformation(iorder);
    end
    transformation=type5_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,5,iorder)=transformation(iorder);
    end
    transformation=type6_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,6,iorder)=transformation(iorder);
    end
    transformation=type7_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,7,iorder)=transformation(iorder);
    end
    transformation=type8_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation_list(iorientation,:));
    for iorder=1:1:7
        sol(iorientation,8,iorder)=transformation(iorder);
    end
end

end

function sol=type1_transformer_specific_3d(alpha,beta,gamma,eX,eY,eZ,orientation)
sol=zeros(1,7);
sol(1)=...
    local_basis_3d(alpha,beta,gamma,...
    orientation,[0 0 0],[0 0 0]);
df=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[0 0 0],eX,eY,eZ,1:1:6);
for idf=1:1:6
    sol(idf+1)=df(idf);
end
end

function sol=type2_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus=alpha*gps;
summation1=0;
summation2=0;
summation3=0;
summation4=0;
summation5=0;
for ialpha=1:1:3
    summation1=summation1+w(ialpha)...
        *local_basis_3d(gps_plus(ialpha),beta,gamma,...
        orientation,[1 0 0],[1 0 0]);
    df=get_fg_dervs_3d(gps_plus(ialpha),beta,gamma,...
        orientation,[1 0 0],eX,eY,eZ,[2 3 5 6]);
    summation2=summation2+w(ialpha)*df(1);
    summation3=summation3+w(ialpha)*df(2);
    summation4=summation4+w(ialpha)*df(3);
    summation5=summation5+w(ialpha)*df(4);
    
    %------------------------------------------
%     if orientation(1)==1 && orientation(2)==0 && orientation(3)==0
%         [df(1) w(ialpha)*df(1)]
%         
%         gps_plus(ialpha)
%         df(1)
%         
%         error('warui')
%     end
%     par=10;
%     for ijk=1
%         pt1=(1/par)/2;
%         pt2=1-pt1;
%         pts=linspace(pt1,pt2,par);
%         summation_=0;
%         for ijk2=1:1:par
%             df=get_fg_dervs_3d(pts(ijk2),beta,gamma,...
%                 orientation,[1 0 0],eX,eY,eZ,2);
%             summation_=summation_+df;
%         end
%         summation_/par
%         
%         
%     end
    %-------------------------------------------------------
    
end
mult=alpha*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eX,[1 0 0]);
sol(1)=summation1*mult;
sol(3)=summation2*mult;

%--------------------------------------
% if orientation(1)==1 && orientation(2)==0 && orientation(3)==0
%     sol(3)
%     mult
%     
%     par=1000;
%     pt1=(1/par)/2;
%     pt2=1-pt1;
%     dist=1/par;
%     pts=linspace(pt1,pt2,par);
%     summation_=0;
%     for ijk2=1:1:par
%         df=get_fg_dervs_3d(pts(ijk2),0,0,...
%             orientation,[1 0 0],eX,eY,eZ,2);
%         summation_=summation_+df*dist;
%     end
%     summation_
%         
%         summation2
%     
%     error('Watch here')
% end
%----------------------------------------------

sol(4)=summation3*mult;
sol(6)=summation4*mult;
sol(7)=summation5*mult;
sol(2)=local_basis_3d(alpha,beta,gamma,...
    orientation,[1 0 0],[1 0 0]);
sol(5)=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[1 0 0],eX,eY,eZ,1);
end

function sol=type3_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus=beta*gps;
summation1=0;
summation2=0;
summation3=0;
summation4=0;
summation5=0;
for ibeta=1:1:3
    summation1=summation1+w(ibeta)...
        *local_basis_3d(alpha,gps_plus(ibeta),gamma,...
        orientation,[0 1 0],[0 1 0]);
    df=get_fg_dervs_3d(alpha,gps_plus(ibeta),gamma,...
        orientation,[0 1 0],eX,eY,eZ,[1 3 4 6]);
    summation2=summation2+w(ibeta)*df(1);
    summation3=summation3+w(ibeta)*df(2);
    summation4=summation4+w(ibeta)*df(3);
    summation5=summation5+w(ibeta)*df(4);
end
mult=beta*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eY,[0 1 0]);
sol(1)=summation1*mult;
sol(2)=summation2*mult;
sol(4)=summation3*mult;
sol(5)=summation4*mult;
sol(7)=summation5*mult;
sol(3)=local_basis_3d(alpha,beta,gamma,...
    orientation,[0 1 0],[0 1 0]);
sol(6)=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[0 1 0],eX,eY,eZ,2);
end

function sol=type4_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus=gamma*gps;
summation1=0;
summation2=0;
summation3=0;
summation4=0;
summation5=0;
for igamma=1:1:3
    summation1=summation1+w(igamma)...
        *local_basis_3d(alpha,beta,gps_plus(igamma),...
        orientation,[0 0 1],[0 0 1]);
    df=get_fg_dervs_3d(alpha,beta,gps_plus(igamma),...
        orientation,[0 0 1],eX,eY,eZ,[1 2 4 5]);
    summation2=summation2+w(igamma)*df(1);
    summation3=summation3+w(igamma)*df(2);
    summation4=summation4+w(igamma)*df(3);
    summation5=summation5+w(igamma)*df(4);
end
mult=gamma*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eZ,[0 0 1]);
sol(1)=summation1*mult;
sol(2)=summation2*mult;
sol(3)=summation3*mult;
sol(5)=summation4*mult;
sol(6)=summation5*mult;
sol(4)=local_basis_3d(alpha,beta,gamma,...
    orientation,[0 0 1],[0 0 1]);
sol(7)=get_fg_dervs_3d(alpha,beta,gamma,...
    orientation,[0 0 1],eX,eY,eZ,3);
end

function sol=type5_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_alpha=alpha*gps;
gps_plus_beta=beta*gps;
summation=0;
for ibeta=1:1:3
    summation=summation+w(ibeta)...
        *local_basis_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[1 1 0],[1 1 0]);
end
sol(2)=summation*beta...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
summation=0;
for ialpha=1:1:3
    summation=summation+w(ialpha)...
        *local_basis_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 1 0],[1 1 0]);
end
sol(3)=summation*alpha...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
summation1=0;
summation2=0;
summation3=0;
for ialpha=1:1:3
    for ibeta=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fz=get_fg_dervs_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,...
            orientation,[1 1 0],eX,eY,eZ,3);
        summation1=summation1+w(ialpha)*w(ibeta)*determinant*fz;
        fzz=get_fg_dervs_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,...
            orientation,[1 1 0],eX,eY,eZ,6);
        summation2=summation2+w(ialpha)*w(ibeta)*determinant*fzz;
        summation3=summation3+w(ialpha)*w(ibeta)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,...
            orientation,[1 1 0],[1 1 0]);
    end
end
sol(4)=summation1*alpha*beta;
sol(7)=summation2*alpha*beta;
sol(1)=summation3*alpha*beta;
summation=0;
for ibeta=1:1:3
    fx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[1 1 0],eX,eY,eZ,1);
    summation=summation+w(ibeta)*fx;
end
sol(5)=summation*beta...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eY,[0 1 0]);
summation=0;
for ialpha=1:1:3
    fy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 1 0],eX,eY,eZ,2);
    summation=summation+w(ialpha)*fy;
end
sol(6)=summation*alpha...
    *map_to_global_1_compt_3d(alpha,beta,gamma,eX,[1 0 0]);
end

function sol=type6_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_alpha=alpha*gps;
gps_plus_gamma=gamma*gps;
summation1=0;
summation2=0;
summation3=0;
for ialpha=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[0 0 1]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),...
            orientation,[1 0 1],eX,eY,eZ,2);
        summation1=summation1+w(ialpha)*w(igamma)*determinant*fy;
        fyy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),...
            orientation,[1 0 1],eX,eY,eZ,5);
        summation2=summation2+w(ialpha)*w(igamma)*determinant*fyy;
        summation3=summation3+w(ialpha)*w(igamma)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),...
            orientation,[1 0 1],[1 0 1]);
    end
end
mult=alpha*gamma;
sol(3)=summation1*mult;
sol(6)=summation2*mult;
sol(1)=summation3*mult;
summation1=0;
summation2=0;
for igamma=1:1:3
    summation1=summation1+w(igamma)...
        *local_basis_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[1 0 1],[1 0 1]);
    fx=get_fg_dervs_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[1 0 1],eX,eY,eZ,1);
    summation2=summation2+w(igamma)*fx;
end
mult=gamma*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eZ,[0 0 1]);
sol(2)=summation1*mult;
sol(5)=summation2*mult;
summation1=0;
summation2=0;
for ialpha=1:1:3
    summation1=summation1+w(ialpha)...
        *local_basis_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 0 1],[1 0 1]);
    fz=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gamma,...
        orientation,[1 0 1],eX,eY,eZ,3);
    summation2=summation2+w(ialpha)*fz;
end
mult=alpha*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eX,[1 0 0]);
sol(4)=summation1*mult;
sol(7)=summation2*mult;
end

function sol=type7_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_beta=beta*gps;
gps_plus_gamma=gamma*gps;
summation1=0;
summation2=0;
summation3=0;
for ibeta=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
        term2=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
        term4=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
        determinant=term1*term2-term3*term4;
        fx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
            orientation,[0 1 1],eX,eY,eZ,1);
        summation1=summation1+w(ibeta)*w(igamma)*determinant*fx;
        fxx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
            orientation,[0 1 1],eX,eY,eZ,4);
        summation2=summation2+w(ibeta)*w(igamma)*determinant*fxx;
        summation3=summation3+w(ibeta)*w(igamma)*determinant...
            *local_basis_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
            orientation,[0 1 1],[0 1 1]);
    end
end
mult=beta*gamma;
sol(2)=summation1*mult;
sol(5)=summation2*mult;
sol(1)=summation3*mult;
summation1=0;
summation2=0;
for igamma=1:1:3
    summation1=summation1+w(igamma)...
        *local_basis_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[0 1 1],[0 1 1]);
    fy=get_fg_dervs_3d(alpha,beta,gps_plus_gamma(igamma),...
        orientation,[0 1 1],eX,eY,eZ,2);
    summation2=summation2+w(igamma)*fy;
end
mult=gamma*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eZ,[0 0 1]);
sol(3)=summation1*mult;
sol(6)=summation2*mult;
summation1=0;
summation2=0;
for ibeta=1:1:3
    summation1=summation1+w(ibeta)...
        *local_basis_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[0 1 1],[0 1 1]);
    fz=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gamma,...
        orientation,[0 1 1],eX,eY,eZ,3);
    summation2=summation2+w(ibeta)*fz;
end
mult=beta*map_to_global_1_compt_3d(...
    alpha,beta,gamma,eY,[0 1 0]);
sol(4)=summation1*mult;
sol(7)=summation2*mult;
end

function sol=type8_transformer_specific_3d(alpha,beta,gamma,gps,w,eX,eY,eZ,orientation)
sol=zeros(1,7);
gps_plus_alpha=alpha*gps;
gps_plus_beta=beta*gps;
gps_plus_gamma=gamma*gps;
summation=0;
J=zeros(3,3);
for ialpha=1:1:3
    for ibeta=1:1:3
        for igamma=1:1:3
            J(1,1)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eX,[1 0 0]);
            J(1,2)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eX,[0 1 0]);
            J(1,3)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eX,[0 0 1]);
            J(2,1)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[1 0 0]);
            J(2,2)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
            J(2,3)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
            J(3,1)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[1 0 0]);
            J(3,2)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
            J(3,3)=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
            summation=summation+w(ialpha)*w(ibeta)*w(igamma)*det(J)...
                *local_basis_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gps_plus_gamma(igamma),...
                orientation,[1 1 1],[1 1 1]);
        end
    end
end
sol(1)=summation*alpha*beta*gamma;
summation=0;
for ibeta=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
        term2=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
        term4=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
        determinant=term1*term2-term3*term4;
        summation=summation+w(ibeta)*w(igamma)*determinant...
            *local_basis_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),orientation,[1 1 1],[1 1 1]);
    end
end
sol(2)=summation*beta*gamma;
summation=0;
for ialpha=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[0 0 1]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[1 0 0]);
        determinant=term1*term2-term3*term4;
        summation=summation+w(ialpha)*w(igamma)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),orientation,[1 1 1],[1 1 1]);
    end
end
sol(3)=summation*alpha*gamma;
summation=0;
for ialpha=1:1:3
    for ibeta=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[1 0 0]);
        determinant=term1*term2-term3*term4;
        summation=summation+w(ialpha)*w(ibeta)*determinant...
            *local_basis_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,orientation,[1 1 1],[1 1 1]);
    end
end
sol(4)=summation*alpha*beta;
summation=0;
for ibeta=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 1 0]);
        term2=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eY,[0 0 1]);
        term4=map_to_global_1_compt_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),eZ,[0 1 0]);
        determinant=term1*term2-term3*term4;
        fx=get_fg_dervs_3d(alpha,gps_plus_beta(ibeta),gps_plus_gamma(igamma),orientation,[1 1 1],eX,eY,eZ,1);
        summation=summation+w(ibeta)*w(igamma)*determinant*fx;
    end
end
sol(5)=summation*beta*gamma;
summation=0;
for ialpha=1:1:3
    for igamma=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[0 0 1]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eX,[0 0 1]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),eZ,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fy=get_fg_dervs_3d(gps_plus_alpha(ialpha),beta,gps_plus_gamma(igamma),orientation,[1 1 1],eX,eY,eZ,2);
        summation=summation+w(ialpha)*w(igamma)*determinant*fy;
    end
end
sol(6)=summation*alpha*gamma;
summation=0;
for ialpha=1:1:3
    for ibeta=1:1:3
        term1=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[1 0 0]);
        term2=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[0 1 0]);
        term3=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eX,[0 1 0]);
        term4=map_to_global_1_compt_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,eY,[1 0 0]);
        determinant=term1*term2-term3*term4;
        fz=get_fg_dervs_3d(gps_plus_alpha(ialpha),gps_plus_beta(ibeta),gamma,orientation,[1 1 1],eX,eY,eZ,3);
        summation=summation+w(ialpha)*w(ibeta)*determinant*fz;
    end
end
sol(7)=summation*alpha*beta;
end