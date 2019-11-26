function [logS_export,limS_export,E_updated_export,...
    max_c_export,min_c_export,intensity_specific]=plots_output(...
    fig_output_dim,dpi,axis_font,font_size,...
    output_path,export_fig,...
    c,nx,ny,alpha,beta,fc,ci_ave,...
    iframe,t,logS_export,limS_export,xlen,ylen,limSstretcher,...
    xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
    frequency_domain,...
    coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,...
    E_updated_export,max_c_export,min_c_export,diff,...
    export_spec,k_specific,intensity_specific,...
    nworkers)


if export_spec==1
    
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    [mag_spec,~]=fourier_analysis_2d(alpha,beta,fc,c_nodal,ci_ave,0,nworkers);
    maxS=max(mag_spec);

    total_energy=evaluate_total_energy(...
        coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
        ny,c,weights,diff,...
        grad_T);

    if iframe==1
        E_updated_export=total_energy;
        max_c_export=max(max(c_nodal));
        min_c_export=min(min(c_nodal));
        
    else
        E_updated_export(iframe)=total_energy;
        max_c_export(iframe)=max(max(c_nodal));
        min_c_export(iframe)=min(min(c_nodal));
    end


    if iframe==1
        logS_export=[t log(maxS)];
    else
        logS_export(iframe,:)=[t log(maxS)];
    end
    while limS_export<=maxS
        limS_export=limS_export*limSstretcher;
    end


    fig=figure(export_fig);
    % set(fig,'Position',fig_output_dim)
    set(fig,'PaperUnits','inches','PaperPosition',fig_output_dim)
    % axes('FontSize',30,'FontName',axis_font)

    x_coord=linspace(0,xlen,nx);
    y_coord=linspace(0,ylen,ny);

    % figure(1988)
    % tru=reshape(c_nodal,[1,nx*ny]);
    % pd_domain=linspace(0,1,100);
    % pd=pdf(makedist('Normal'),pd_domain);
    % % plot(pd_domain,pd)
    % histogram(tru,100)
    % xlim([0 1])
    % print('-dpng', strcat(output_path,'figures/profile/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))
    % return
    % error('adfasdf')


    if iframe==1
        %========================================== LEFT OFF
    %     illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
        if length(T)==1
            if length(ci)==1
                h=plot(xspinodal,yspinodal,'k',xbinodal,ybinodal,'k',ci,T,'k*');
            elseif length(ci)==2
                h=plot(xspinodal,yspinodal,'k',xbinodal,ybinodal,'k',ci,[T T],'k*');
            end
        elseif length(T)==2
            if length(ci)==1
                h=plot(xspinodal,yspinodal,'k',xbinodal,ybinodal,'k',[ci ci],T,'k*');
            elseif length(ci)==2
                h=plot(xspinodal,yspinodal,'k',xbinodal,ybinodal,'k',ci,T,'k*');
            end
        end
        set(h(2),'LineWidth',2)
        legend('Spinodal curve','Binodal curve','Point of interest','Location','southeast');
        grid on
        grid minor
        box on
    %     title('Simulation target')
    %     saveas(gcf,[output_path 'figures/simulation_range/phase_diagram.png'])
        xlabel('Temperature, T')
        ylabel('Concentration, c')
        set(gca,'FontName',axis_font,'FontSize',font_size)

        %=============================================
        print('-dpng', strcat(output_path,'figures/simulation_range/phase_diagram.png') , strcat('-r',num2str(dpi)))
    end

    % subplot(1,2,1)

    range_=[min(min(c_nodal)),max(max(c_nodal))];
    % surf(x_coord,y_coord,c_nodal,'edgecolor','none');
    surf(x_coord,y_coord,c_nodal,'edgecolor','none');
    % shading interp
    axis([0 1 0 1 0 1])
    % xlabel('x'),ylabel('y'),zlabel('Concentration, u')
    try
        caxis(range_)
    catch
        %Forget about adjusting axis
    end
    grid on
    grid minor
    box on
    xlabel('x')
    ylabel('y')
    zlabel('Concentration, c')

    camlight
    % subplot(1,2,2)

    % contourf(x_coord,y_coord,c_nodal,[0 ci_ave],'linestyle', 'none');
    % colormap([0 0 0;1 1 1])
    % xlabel('x');ylabel('y')
    % grid on
    % grid minor
    % box on
    % axis equal
    % xlabel('x')
    % ylabel('y')
    % zlabel('Concentration, u')
    % camlight
    % title('Surface plot')
    % saveas(gcf,[output_path 'figures/profile/' num2str(iframe) '.png'])

    % h=colorbar('location','eastoutside');
    % yt=get(h,'XTick');
    % 
    % yt'
    % yt=yt-ci_ave*ones(1,length(yt))
    % ci_ave
    % 
    % yt=[0 1 2]
    % 
    % set(h,'XTickLabel',yt);

    % set(gca,'FontName',axis_font,'FontSize',font_size)

    % error('dfafdafd')
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/profile/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))




    contourf(x_coord,y_coord,c_nodal,[0 ci_ave],'linestyle', 'none');
    colormap([0 0 0;1 1 1])
    xlabel('x');ylabel('y')
    box on
    % saveas(gcf,[output_path 'figures/phases/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/phases/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))


    plot(frequency_domain,mag_spec,'k')
    ytickformat('%,.1f')
    xlabel('Frequency, k')
    ylabel('|| {\fontname{Lucida Calligraphy}F} ||^2')
    ylim([0 limS_export])
    grid on
    grid minor
    box on
    % saveas(gcf,[output_path 'figures/fourier_transform/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/fourier_transform/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))

    plot(logS_export(:,1),logS_export(:,2),'k*')
    xlim([0 inf])
    grid on
    grid minor
    ytickformat('%,1.1f')
    xlabel('Time, t')
    ylabel('log(|| {{\fontname{Lucida Calligraphy}F_{c}} ||_{max}}^2)')
    box on
    % saveas(gcf,[output_path 'figures/stage_view/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/stage_view/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))

    plot(logS_export(:,1),max_c_export,'k*')
    axis([0 inf 0 1])
    ytickformat('%,1.1f')
    grid on
    grid minor
    % ytickformat('%,1.1f')
    xlabel('Time, t')
    ylabel('Maximum concentration, c_m_a_x')
    box on
    % saveas(gcf,[output_path 'figures/max_c/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/max_c/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))


    plot(logS_export(:,1),min_c_export,'k*')
    axis([0 inf 0 1])
    grid on
    grid minor
    % ytickformat('%,1.1f')
    xlabel('Time, t')
    ylabel('Minimum concentration, c_m_i_n')
    box on
    % saveas(gcf,[output_path 'figures/min_c/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/min_c/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))

    %----------

    plot(logS_export(:,1),E_updated_export-ones(1,iframe)*E_updated_export(1),'k*')
    
    % xlim([0 inf])
    axis([0 inf -inf 0])
    xlabel('Time, t');ylabel('Energy change, \DeltaE')
    grid on
    grid minor
    ytickformat('%,.1f')
    box on

    % saveas(gcf,[output_path 'figures/energy/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/energy/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))


    close(fig)



    % mkdir([output_path 'profile/'])
    %     mkdir([output_path 'phases/'])
    %     mkdir([output_path 'fourier_transform/'])
    %     mkdir([output_path 'stage_view/'])
    %     mkdir([output_path 'energy/'])
    
elseif export_spec==2
    
    fig=figure(export_fig);
    potential=evaluate_potential(c,nex,ney,ny,n1,n2,Two_chi_n1,...
        weights,grad_T,coef_T);
    dx_2=xlen/(2*nex);
    dy_2=ylen/(2*ney);
    x_coords=linspace(dx_2,xlen-dx_2,nex);
    y_coords=linspace(dy_2,ylen-dy_2,ney);
    surf(x_coords,y_coords,potential,'edgecolor','none')
    xlabel('x');ylabel('y');zlabel('\mu_2-\mu_1')
    xlim([0 xlen])
    ylim([0 ylen])
    camlight
    grid on
    grid minor
    ytickformat('%,.1f')
    box on
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,'figures/potential/',num2str(iframe),'.png') , strcat('-r',num2str(dpi)))
    close(fig)
    
elseif export_spec==3
    
    fig=figure(export_fig);
    %----------------
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    [mag_spec,~]=fourier_analysis_2d(2,beta,k_specific,c_nodal,ci_ave,0,nworkers);
    intensity=log(mag_spec(2));
    %------------------
    if iframe>1
        intensity_specific(iframe,:)=[t intensity];
    elseif iframe==1
        intensity_specific=[0 intensity];
    end
    plot(intensity_specific(:,1),intensity_specific(:,2),'k*')
    xlim([0 inf])
    grid on
    grid minor
    ytickformat('%,1.1f')
    xlabel('Time, t')
%     sss=num2str(k_specific);
    ylabel(['Log of squared magnitude, log(|| {\itF}' '_k_=_' '{' num2str(k_specific) '}' ' ||^2)'])
    box on
    % saveas(gcf,[output_path 'figures/stage_view/' num2str(iframe) '.png'])
    set(gca,'FontName',axis_font,'FontSize',font_size)
    print('-dpng', strcat(output_path,['figures/k=' num2str(k_specific) '/'],num2str(iframe),'.png') , strcat('-r',num2str(dpi)))
    close(fig)
    
end


end