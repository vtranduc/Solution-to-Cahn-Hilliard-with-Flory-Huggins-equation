function [frame_cap,E_updated,ac_t_updated,logS,limS]=visual_2d(...
    figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,t,...
    frequency_domain,ci_ave,limS,logS,...
    xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
    iframe,...
    ac_t,...
    alpha,beta,fc,limSstretcher,...
    nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size,...
    nworkers)

% function [frame_cap,E_updated,ac_t_updated,logS,limS]=visual_2d_2018(...
%     figure_analysis,c,nx,ny,x_coord,y_coord,minz,maxz,t,...
%     frequency_domain,ci_ave,limS,logS,...
%     xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
%     ne,nex,ney,weights,coef_T,Two_chi_n1,n1,n2,...
%     dxdy,iframe,E,...
%     co,dt,ac_t,...
%     alpha,beta,fc,limSstretcher,...
%     nx_grad,ny_grad,x_grad,y_grad,contour_lgd_font_size,...
%     grad_T,diff,entropy)

%======================ANALYSIS============================================

if figure_analysis==1
    
    error('Energy and potential are wrong for figure_analysis==1')
    
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    [DX,DY]=extract_gradient_2D(c,nx,ny);
    laplacian=approximate_laplacian_2d(DX,DY,nx,ny,nex,ney);
    
    %POTENTIAL DOES NOT TAKE INTO ACCOUNT!!!===============================
    [energy,potential]=compute_energy_potential(c_nodal,coef_T,Two_chi_n1,nx,ny,n1,n2,DX,DY,grad_T,diff,T,entropy);
    %================================================================
    
    
    [mag_spec,phase_spec]=fourier_analysis_2d(alpha,beta,fc,c_nodal,ci_ave,1);
    maxS=max(mag_spec);

    if t~=0
        logS(iframe,:)=[t log(maxS)];
    elseif t==0
        logS=[t log(maxS)];
    end
    while limS<=maxS
        limS=limS*limSstretcher;
    end
    E_updated=E;
    %=============================================
%     E_updated(iframe)=...
%         compute_total_energy_2d(...
%         coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%         ny,c,weights,...
%         grad_T,diff);
    E_updated(iframe)=...
        evaluate_total_energy(...
        coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
        ny,c,weights,diff,...
        grad_T);
    %===============================================
    growth=(c_nodal-extract_nodal_weights_2D(co,nx,ny))/dt;
    allen_cahn=laplacian-potential;
    
elseif figure_analysis==2
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
elseif figure_analysis==3
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    [mag_spec,~]=fourier_analysis_2d(alpha,beta,fc,c_nodal,ci_ave,0,nworkers);
    maxS=max(mag_spec);
    
    while limS<=maxS
        limS=limS*limSstretcher;
    end
    
    %================================ TEST
%     [maxS,~]=fourier_analysis_2d(2,beta,5.86615,c_nodal,ci_ave,1);
%     
%     maxS=maxS(2);
%     
%     warning('fdafasd')
    %=========================================
    
    if t~=0
        logS(iframe,:)=[t log(maxS)];
    elseif t==0
        logS=[t log(maxS)];
    end
    
    
    
    
end

if figure_analysis==1
    
    %======================MAIN=============================================

    h1=subplot(3,6,[1 2 7 8]);
    h1_pos=get(h1,'Position');
    main_profile(c_nodal,x_coord,y_coord,minz,maxz)

    %=======================================================================

    %=====================PHASE DIAGRAM=======================================

    subplot(3,6,13)
    phase_profile(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)

    %========================================================================

    %=====================TWO PHASES=======================================

    subplot(3,6,14)
    two_phases_profile(x_coord,y_coord,c_nodal,ci_ave)
    
    %========================================================================

    %=====================Growth=======================================
    
    subplot(3,6,3)
    growth_profile(t,x_coord,y_coord,growth,contour_lgd_font_size)

    %========================================================================

    %======================ALLEN CAHN==========================================

    subplot(3,6,9)
    allen_cahn_profile(x_coord,y_coord,allen_cahn,contour_lgd_font_size)
    
    subplot(3,6,15)
    ac_t_updated=min_max_ave_allen_cahn_profile(ac_t,t,iframe,allen_cahn,logS);

    %========================================================================

    %======================CONCENTRATION=======================================

    subplot(3,6,4);
    concentration_profile(x_coord,y_coord,c_nodal,contour_lgd_font_size)

    subplot(3,6,10)
    gradient_profile(nx_grad,ny_grad,DX,DY,x_coord,y_coord,c_nodal,...
        x_grad,y_grad,contour_lgd_font_size)
    
    subplot(3,6,16)
    laplacian_profile(t,x_coord,y_coord,laplacian,contour_lgd_font_size)

    %========================================================================

    %======================FOURIER ANALYSIS====================================
    
    subplot(3,6,5)
    fourier_magnitude_profile(frequency_domain,mag_spec,limS)
    
    subplot(3,6,11)
    fourier_phase_profile(frequency_domain,phase_spec)
    
    subplot(3,6,17)
    fourier_logged_magnitude_profile(logS)

    %========================================================================

    %======================ENERGY=============================================

    subplot(3,6,6)
    energy_profile(x_coord,y_coord,energy,contour_lgd_font_size)
    
    subplot(3,6,12)
    potential_profile(x_coord,y_coord,potential,contour_lgd_font_size)

    subplot(3,6,18)
    total_energy_profile(logS,E_updated)

    %=====================================================================

    universal_colorbar(h1_pos,[0.2 0.9])

%======================================================================
%======================================================================
%======================================================================
%======================================================================
%======================================================================
%==========DIFFERENT ANALYSIS==========================================
%======================================================================
%======================================================================
%======================================================================
%======================================================================
%======================================================================
    
elseif figure_analysis==2
    main_profile(c_nodal,x_coord,y_coord,minz,maxz)
    colorbar
    E_updated=NaN;ac_t_updated=NaN;logS=NaN;limS=NaN;
elseif figure_analysis==3
    subplot(2,4,[1,2,5,6])
    main_profile(c_nodal,x_coord,y_coord,minz,maxz)
    colorbar('westoutside');
    subplot(2,4,3)
    phase_profile(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
    grid on
    grid minor
    subplot(2,4,4)
    two_phases_profile(x_coord,y_coord,c_nodal,ci_ave)
    subplot(2,4,7)
    fourier_magnitude_profile(frequency_domain,mag_spec,limS)
    subplot(2,4,8)
    fourier_logged_magnitude_profile(logS)
    E_updated=NaN;ac_t_updated=NaN;
end

str=sprintf('Simulation of Cahn Hilliard in 2D\n t=%d',t);
suptitle(str)

frame_cap=getframe(gcf);

end

%======================================================================
%======================================================================
%======================================================================
%======================================================================
%======================================================================
%==========PLOTTING FUNCTIONS==========================================
%======================================================================
%======================================================================
%======================================================================
%======================================================================
%======================================================================

function main_profile(c_nodal,x_coord,y_coord,minz,maxz)
range_=[min(min(c_nodal)),max(max(c_nodal))];
surf(x_coord,y_coord,c_nodal,'edgecolor','none');
% shading interp
axis([0 1 0 1 minz maxz])
xlabel('x'),ylabel('y'),zlabel('c')
try
    caxis(range_)
catch
    %Forget about adjusting axis
end
camlight
title('Surface plot')
end

function phase_profile(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
grid on
grid minor
title('Simulation target')
end

function growth_profile(t,x_coord,y_coord,growth,contour_lgd_font_size)
if t~=0
    contourf(x_coord,y_coord,growth,'edgecolor','none');
elseif t==0
    warning('off','all')
    contourf(x_coord,y_coord,growth,'edgecolor','none');
    warning('on','all')
end
xlabel('x'),ylabel('y')
title('\partialc/\partialt')
str=sprintf('Max=%.1e,Min=%.1e',max(max(growth)),min(min(growth)));
l=legend(str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function allen_cahn_profile(x_coord,y_coord,allen_cahn,contour_lgd_font_size)
contourf(x_coord,y_coord,allen_cahn,'edgecolor','none')
xlabel('x'),ylabel('y')
title('\nabla^2c - \partialf/\partialc')
str=sprintf('Max=%.1e,Min=%.1e',max(max(allen_cahn)),min(min(allen_cahn)));
l=legend(str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function ac_t_updated=min_max_ave_allen_cahn_profile(ac_t,t,iframe,allen_cahn,logS)
ac_t_updated=ac_t;
ac_t_updated(iframe,:)=[max(max(allen_cahn)) min(min(allen_cahn)) mean(mean(allen_cahn))];
if t==0
    plot(logS(:,1),ac_t_updated(:,1),logS(:,1),ac_t_updated(:,2),logS(:,1),ac_t_updated(:,3))
else
    [nsteps,~]=size(logS);
    plot(logS(2:nsteps,1),ac_t_updated(2:nsteps,1),logS(2:nsteps,1),ac_t_updated(2:nsteps,2),logS(2:nsteps,1),ac_t_updated(2:nsteps,3))
end
grid on
grid minor
ytickformat('%,.1f')
xlabel('Time, t')
ylabel('Allen-Cahn value')
title('Max, min, mean of \nabla^2c - \partialf/\partialc')
end

function two_phases_profile(x_coord,y_coord,c_nodal,ci_ave)
contourf(x_coord,y_coord,c_nodal,[0 ci_ave],'linestyle', 'none');
xlabel('x');ylabel('y')
title({'Two phases as distinguised','by average concentration'})
end

function concentration_profile(x_coord,y_coord,c_nodal,contour_lgd_font_size)
contourf(x_coord,y_coord,c_nodal,'edgecolor','none');
axis([0 1 0 1])
xlabel('x');ylabel('y')
title('c')
str=sprintf('Max=%.1e,Min=%.1e',max(max(c_nodal)),min(min(c_nodal)));
l=legend(str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function gradient_profile(nx_grad,ny_grad,DX,DY,x_coord,y_coord,c_nodal,...
    x_grad,y_grad,contour_lgd_font_size)
nnx=length(nx_grad);
nny=length(ny_grad);
DX_=zeros(nny,nnx);
DY_=zeros(nny,nnx);
for inx=1:1:nnx
    for iny=1:1:nny
        DX_(iny,inx)=DX(ny_grad(iny),nx_grad(inx));
        DY_(iny,inx)=DY(ny_grad(iny),nx_grad(inx));
    end
end
contour(x_coord,y_coord,c_nodal)
hold on
h=quiver(x_grad,y_grad,DX_,DY_,'color','black');
hold off
grid on
grid minor
xlabel('x');ylabel('y')
title('\nablac')
abs_len=(DX_.^2+DY_.^2).^0.5;
str=sprintf('Max=%.1e,Min=%.1e',max(max(abs_len)),min(min(abs_len)));
l=legend(h,str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function laplacian_profile(t,x_coord,y_coord,laplacian,contour_lgd_font_size)
if t~=0
    contourf(x_coord,y_coord,laplacian,'edgecolor','none')
else
    warning('off','all')
    contourf(x_coord,y_coord,laplacian,'edgecolor','none')
    warning('on','all')
end
xlabel('x');ylabel('y')
title('\nabla^2c')
str=sprintf('Max=%.1e,Min=%.1e',max(max(laplacian)),min(min(laplacian)));
l=legend(str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function fourier_magnitude_profile(frequency_domain,mag_spec,limS)
plot(frequency_domain,mag_spec)
ytickformat('%,.1f')
xlabel('Frequency domain, k')
ylabel('||\it F ||^2')
ylim([0 limS])
grid on
grid minor
title('||\it F ||^2')
end

function fourier_phase_profile(frequency_domain,phase_spec)
plot(frequency_domain,phase_spec)
ytickformat('%,.1f')
xlabel('Frequency domain, k')
ylabel('Phase')
grid on
grid minor
title('Phase')
end

function fourier_logged_magnitude_profile(logS)
plot(logS(:,1),logS(:,2),'*')
grid on
grid minor
ytickformat('%,1.1f')
xlabel('Time, t')
ylabel('log(|| {\itF} ||^2)')
title('log(|| {\itF} ||^2) vs Time')
end

function energy_profile(x_coord,y_coord,energy,contour_lgd_font_size)
contourf(x_coord,y_coord,energy,'edgecolor','none')
xlabel('x'),ylabel('y')
title('Energy, E_v*L^2/\kappa')
str=sprintf('Max=%.1e,Min=%.1e',max(max(energy)),min(min(energy)));
l=legend(str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function potential_profile(x_coord,y_coord,potential,contour_lgd_font_size)
contourf(x_coord,y_coord,potential,'edgecolor','none')
xlabel('x'),ylabel('y')
title('Potential, P_v*L^2/\kappa')
str=sprintf('Max=%.1e,Min=%.1e',max(max(potential)),min(min(potential)));
l=legend(str,'Location','southwest');
l.FontSize=contour_lgd_font_size;
end

function total_energy_profile(logS,E_updated)
plot(logS(:,1),E_updated,'*')
xlabel('Time');ylabel('Total Energy, E*L^2/\kappa')
grid on
grid minor
ytickformat('%,.1f')
title('Energy vs time')
end

function universal_colorbar(h1_pos,ratio)
bar_size=[h1_pos(1) 1].*ratio;
colorbar('Position', [(h1_pos(1)-bar_size(1))/2 0.5-bar_size(2)/2 bar_size(1) bar_size(2)],'YTick',[])
end

function [E,P]=compute_energy_potential(c_nodal,coef_T,Two_chi_n1,nx,ny,n1,n2,DX,DY,grad_T,diff,T,entropy)
E=zeros(ny,nx);
P=zeros(ny,nx);
chi_n1=Two_chi_n1/2;
if grad_T==0
    two_diffT=2*coef_T;
    for inx=1:1:nx
        for iny=1:1:ny
            c2=1-c_nodal(iny,inx);
            c1logged=log(c_nodal(iny,inx));
            c2logged=log(c2);
            E(iny,inx)=two_diffT*(c_nodal(iny,inx)*c1logged/n1+c2*c2logged/n2+c_nodal(iny,inx)*c2*chi_n1)+DX(iny,inx)^2+DY(iny,inx)^2;
            P(iny,inx)=two_diffT*((c1logged+1)/n1-(c2logged+1)/n2+chi_n1*(1-2*c_nodal(iny,inx)));
        end
    end
elseif grad_T==1
    two_diffT=linspace(T(1),T(2),nx);
    chi_n1=get_Two_chi_n1_with_order(two_diffT,entropy,n1,0)/2;
    two_diffT=two_diffT*diff*2;
    for inx=1:1:nx
        for iny=1:1:ny
            c2=1-c_nodal(iny,inx);
            c1logged=log(c_nodal(iny,inx));
            c2logged=log(c2);
            E(iny,inx)=two_diffT(inx)*(c_nodal(iny,inx)*c1logged/n1+c2*c2logged/n2+c_nodal(iny,inx)*c2*chi_n1(inx))+DX(iny,inx)^2+DY(iny,inx)^2;
            P(iny,inx)=two_diffT(inx)*((c1logged+1)/n1-(c2logged+1)/n2+chi_n1(inx)*(1-2*c_nodal(iny,inx)));
        end
    end
end
end