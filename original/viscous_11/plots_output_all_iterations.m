function plots_output_all_iterations(...
    fig_output_dim,dpi,axis_font,font_size,...
    simulation_folder,export_fig,...
    nx,ny,alpha,beta,fc,ci_ave,...
    xlen,ylen,limSstretcher,...
    xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
    frequency_domain,limS,...
    coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,...
    nframes_last_sim,...
    export_spec,diff,k_specific,...
    nworkers)

if exist([simulation_folder 'figures/'],'dir')
    if exist([simulation_folder 'old_figures/'],'dir')
        rmdir([simulation_folder 'old_figures/'],'s')
    end
%     movefile([simulation_folder 'figures/'],[simulation_folder 'old_figures/'])
    
    if export_spec==1
        if exist([simulation_folder 'figures/simulation_range/'],'dir')
            movefile([simulation_folder 'figures/simulation_range/'],[simulation_folder 'old_figures/simulation_range/'])
        end
        if exist([simulation_folder 'figures/profile/'],'dir')
            movefile([simulation_folder 'figures/profile/'],[simulation_folder 'old_figures/profile/'])
        end
        if exist([simulation_folder 'figures/phases/'],'dir')
            movefile([simulation_folder 'figures/phases/'],[simulation_folder 'old_figures/phases/'])
        end
        if exist([simulation_folder 'figures/fourier_transform/'],'dir')
            movefile([simulation_folder 'figures/fourier_transform/'],[simulation_folder 'old_figures/fourier_transform/'])
        end
        if exist([simulation_folder 'figures/stage_view/'],'dir')
            movefile([simulation_folder 'figures/stage_view/'],[simulation_folder 'old_figures/stage_view/'])
        end
        if exist([simulation_folder 'figures/energy/'],'dir')
            movefile([simulation_folder 'figures/energy/'],[simulation_folder 'old_figures/energy/'])
        end
        if exist([simulation_folder 'figures/max_c/'],'dir')
            movefile([simulation_folder 'figures/max_c/'],[simulation_folder 'old_figures/max_c/'])
        end
        if exist([simulation_folder 'figures/min_c/'],'dir')
            movefile([simulation_folder 'figures/min_c/'],[simulation_folder 'old_figures/min_c/'])
        end 
            
    elseif export_spec==2
        if exist([simulation_folder 'figures/potential/'],'dir')
            movefile([simulation_folder 'figures/potential/'],[simulation_folder 'old_figures/potential/'])
        end 
    elseif export_spec==3
%         mkdir([simulation_folder 'old_figures'])
        if exist([simulation_folder 'figures/k=' num2str(k_specific) '/'],'dir')
            movefile([simulation_folder 'figures/k=' num2str(k_specific) '/'],[simulation_folder 'old_figures/k=' num2str(k_specific) '/'])
        end
    end
end

try
% if 1==1

    if ~exist([simulation_folder 'figures/'],'dir')
        mkdir([simulation_folder 'figures/'])
    end
    
    if export_spec==1
    
        mkdir([simulation_folder 'figures/simulation_range/'])
        mkdir([simulation_folder 'figures/profile/'])
        mkdir([simulation_folder 'figures/phases/'])
        mkdir([simulation_folder 'figures/fourier_transform/'])
        mkdir([simulation_folder 'figures/stage_view/'])
        mkdir([simulation_folder 'figures/energy/'])
        mkdir([simulation_folder 'figures/max_c/'])
        mkdir([simulation_folder 'figures/min_c/'])
    elseif export_spec==2
        mkdir([simulation_folder 'figures/potential/'])
    elseif export_spec==3
        mkdir([simulation_folder 'figures/k=' num2str(k_specific) '/'])
    end


    output_conc=[simulation_folder 'concentration/iteration_'];
    output_time=[simulation_folder 'time/iteration_'];
    c=dlmread([output_conc '1']);
    time=0;
    [logS_export,limS_export,E_updated_export,...
        max_c_export,min_c_export,intensity_specific]=plots_output(...
        fig_output_dim,dpi,axis_font,font_size,...
        simulation_folder,export_fig,...
        c,nx,ny,alpha,beta,fc,ci_ave,...
        1,time,NaN,limS,xlen,ylen,limSstretcher,...
        xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
        frequency_domain,...
        coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,...
        NaN,NaN,NaN,diff,...
        export_spec,k_specific,NaN,...
        nworkers);

    for iframe=2:1:nframes_last_sim
        c=dlmread([output_conc num2str(iframe)]);
        time=dlmread([output_time num2str(iframe)]);
        [logS_export,limS_export,E_updated_export,...
            max_c_export,min_c_export,intensity_specific]=plots_output(...
            fig_output_dim,dpi,axis_font,font_size,...
            simulation_folder,export_fig,...
            c,nx,ny,alpha,beta,fc,ci_ave,...
            iframe,time,logS_export,limS_export,xlen,ylen,limSstretcher,...
            NaN,NaN,NaN,NaN,NaN,NaN,...
            frequency_domain,...
            coef_T,n1,n2,ne,nex,ney,Two_chi_n1,dxdy,weights,grad_T,...
            E_updated_export,max_c_export,min_c_export,diff,...
            export_spec,k_specific,intensity_specific,...
            nworkers);
        
    end
catch
% elseif 1==2
    warning('Failed to output figures. Old files will be returned if existed')
%     if exist([simulation_folder 'figures/'],'dir')
%         rmdir([simulation_folder 'figures/'],'s')
%     end
%     if exist([simulation_folder 'old_figures/'],'dir')
%         movefile([simulation_folder 'old_figures/'],[simulation_folder 'figures/'])
%     end

    %========================================================


    if export_spec==1
        if exist([simulation_folder 'figures/simulation_range/'],'dir')
            rmdir([simulation_folder 'figures/simulation_range/'],'dir')
        end
        if exist([simulation_folder 'old_figures/simulation_range/'],'dir')
            movefile([simulation_folder 'old_figures/simulation_range/'],[simulation_folder 'figures/simulation_range/'])
        end
        
        if exist([simulation_folder 'figures/profile/'],'dir')
            rmdir([simulation_folder 'figures/profile/'],'dir')
        end
        if exist([simulation_folder 'old_figures/profile/'],'dir')
            movefile([simulation_folder 'old_figures/profile/'],[simulation_folder 'figures/profile/'])
        end
        
        if exist([simulation_folder 'figures/phases/'],'dir')
            rmdir([simulation_folder 'figures/phases/'],'dir')
        end
        if exist([simulation_folder 'old_figures/phases/'],'dir')
            movefile([simulation_folder 'old_figures/phases/'],[simulation_folder 'figures/phases/'])
        end
        
        if exist([simulation_folder 'figures/fourier_transform/'],'dir')
            rmdir([simulation_folder 'figures/fourier_transform/'],'dir')
        end
        if exist([simulation_folder 'old_figures/fourier_transform/'],'dir')
            movefile([simulation_folder 'old_figures/fourier_transform/'],[simulation_folder 'figures/fourier_transform/'])
        end
        
        if exist([simulation_folder 'figures/stage_view/'],'dir')
            rmdir([simulation_folder 'figures/stage_view/'],'dir')
        end
        if exist([simulation_folder 'old_figures/stage_view/'],'dir')
            movefile([simulation_folder 'old_figures/stage_view/'],[simulation_folder 'figures/stage_view/'])
        end
        
        if exist([simulation_folder 'figures/energy/'],'dir')
            rmdir([simulation_folder 'figures/energy/'],'dir')
        end
        if exist([simulation_folder 'old_figures/energy/'],'dir')
            movefile([simulation_folder 'old_figures/energy/'],[simulation_folder 'figures/energy/'])
        end
        
        if exist([simulation_folder 'figures/max_c/'],'dir')
            rmdir([simulation_folder 'figures/max_c/'],'dir')
        end
        if exist([simulation_folder 'old_figures/max_c/'],'dir')
            movefile([simulation_folder 'old_figures/max_c/'],[simulation_folder 'figures/max_c/'])
        end
        
        if exist([simulation_folder 'figures/min_c/'],'dir')
            rmdir([simulation_folder 'figures/min_c/'],'dir')
        end
        if exist([simulation_folder 'old_figures/min_c/'],'dir')
            movefile([simulation_folder 'old_figures/min_c/'],[simulation_folder 'figures/min_c/'])
        end
        
        
%         if exist([simulation_folder 'figures/profile/'],'dir')
%             movefile([simulation_folder 'figures/profile/'],[simulation_folder 'old_figures/profile/'])
%         end
%         if exist([simulation_folder 'figures/phases/'],'dir')
%             movefile([simulation_folder 'figures/phases/'],[simulation_folder 'old_figures/phases/'])
%         end
%         if exist([simulation_folder 'figures/fourier_transform/'],'dir')
%             movefile([simulation_folder 'figures/fourier_transform/'],[simulation_folder 'old_figures/fourier_transform/'])
%         end
%         if exist([simulation_folder 'figures/stage_view/'],'dir')
%             movefile([simulation_folder 'figures/stage_view/'],[simulation_folder 'old_figures/stage_view/'])
%         end
%         if exist([simulation_folder 'figures/energy/'],'dir')
%             movefile([simulation_folder 'figures/energy/'],[simulation_folder 'old_figures/energy/'])
%         end
%         if exist([simulation_folder 'figures/max_c/'],'dir')
%             movefile([simulation_folder 'figures/max_c/'],[simulation_folder 'old_figures/max_c/'])
%         end
%         if exist([simulation_folder 'figures/min_c/'],'dir')
%             movefile([simulation_folder 'figures/min_c/'],[simulation_folder 'old_figures/min_c/'])
%         end 

    elseif export_spec==2
        
        
        if exist([simulation_folder 'figures/potential/'],'dir')
            rmdir([simulation_folder 'figures/potential/'],'dir')
        end
        if exist([simulation_folder 'old_figures/potential/'],'dir')
            movefile([simulation_folder 'old_figures/potential/'],[simulation_folder 'figures/potential/'])
        end
%         if exist([simulation_folder 'figures/potential/'],'dir')
%             movefile([simulation_folder 'figures/potential/'],[simulation_folder 'old_figures/potential/'])
%         end 
    elseif export_spec==3
        
        if exist([simulation_folder 'figures/k=' num2str(k_specific) '/'],'dir')
            rmdir([simulation_folder 'figures/k=' num2str(k_specific) '/'],'dir')
        end
        if exist([simulation_folder 'old_figures/k=' num2str(k_specific) '/'],'dir')
            movefile([simulation_folder 'old_figures/k=' num2str(k_specific) '/'],[simulation_folder 'figures/k=' num2str(k_specific) '/'])
        end
        
%         mkdir([simulation_folder 'old_figures'])
%         if exist([simulation_folder 'figures/k=' num2str(k_specific) '/'],'dir')
%             movefile([simulation_folder 'figures/k=' num2str(k_specific) '/'],[simulation_folder 'old_figures/k=' num2str(k_specific) '/'])
%         end
    end


    %========================================================
    
end

% for iframe=2:1:nframes_last_sim
%     [logS_export,limS_export,E_updated_export,...
%         max_c_export,min_c_export]=plots_output(...
%         fig_output_dim,dpi,axis_font,font_size,...
%         output_path,export_fig,...
%         c,nx,ny,alpha,beta,fc,ci_ave,...
%         iframe,time,logS_export,limS_export,xlen,ylen,limSstretcher,...
%         NaN,NaN,NaN,NaN,NaN,NaN,...
%         frequency_domain,...
%         coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,weights,grad_T,...
%         E_updated_export,max_c_export,min_c_export);
% end


end


% [logS_export,limS_export,E_updated_export,...
%             max_c_export,min_c_export]=plots_output(...
%             fig_output_dim,dpi,axis_font,font_size,...
%             output_path,export_fig,...
%             c,nx,ny,alpha,beta,fc,ci_ave,...
%             iframe,time,NaN,limS,xlen,ylen,limSstretcher,...
%             xspinodal,yspinodal,xbinodal,ybinodal,T,ci,...
%             frequency_domain,...
%             coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,weights,grad_T,...
%             NaN,NaN,NaN);
% 
% 
% [logS_export,limS_export,E_updated_export,...
%                 max_c_export,min_c_export]=plots_output(...
%                 fig_output_dim,dpi,axis_font,font_size,...
%                 output_path,export_fig,...
%                 c,nx,ny,alpha,beta,fc,ci_ave,...
%                 iframe,time,logS_export,limS_export,xlen,ylen,limSstretcher,...
%                 NaN,NaN,NaN,NaN,NaN,NaN,...
%                 frequency_domain,...
%                 coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,weights,grad_T,...
%                 E_updated_export,max_c_export,min_c_export);