function [logS,limS_updated,ac_t,E]=reload_fig_2d(....
    simulation_folder,figure_analysis,iteration_last,nx,ny,...
    alpha,beta,fc,ci_ave,limS,limSstretcher,nworkers)

% function [logS,limS_updated,ac_t,E]=reload_fig_2d(....
%     simulation_folder,figure_analysis,iteration_last,nx,ny,...
%     alpha,beta,fc,ci_ave,limS,limSstretcher,coef_T,n1,n2,ne,...
%     nex,ney,Two_chi_n1,dxdy,weights,...
%     grad_T,diff,T,entropy)

if figure_analysis==0 || figure_analysis==2
    logS=0;E=0;ac_t=[0 0 0];limS_updated=limS;
    return
end
logS=zeros(iteration_last-1,2);
conc_file=[simulation_folder 'concentration/iteration_'];
time_file=[simulation_folder 'time/iteration_'];
if figure_analysis==1
    ac_t=zeros(iteration_last-1,3);
    E=zeros(1,iteration_last-1);
end

f=waitbar(0,'Loading data for figures','Name','Reloading previously obtained data');

for iteration=1:1:iteration_last-1
    
    c=dlmread([conc_file num2str(iteration)]);
    t=dlmread([time_file num2str(iteration)]);
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    
    [mag_spec,~]=fourier_analysis_2d(alpha,beta,fc,c_nodal,ci_ave,0,nworkers);

    logS(iteration,:)=[t log(max(mag_spec))];
    if figure_analysis==1
        error('Currently, feature for figure_analysis==1 is unavailable')
        
%         [DX,DY]=extract_gradient_2D(c,nx,ny);
%         laplacian=approximate_laplacian_2d(DX,DY,nx,ny,nex,ney);
%         potential=compute_potential_2d(...
%             c_nodal,nx,ny,n1,n2,coef_T,diff,grad_T,Two_chi_n1,entropy,T);
%         allen_cahn=laplacian-potential;
%         ac_t(iteration,:)=min_max_ave_allen_cahn(allen_cahn);
%         E(iteration)=...
%             compute_total_energy_2d(...
%             coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%             ny,c,weights,...
%             grad_T,diff);
        
%         total_energy=compute_total_energy_2d(...
%     coef_T,n1,n2,ne,ney,Two_chi_n1,dxdy,...
%     ny,c,weights,...
%     grad_T,diff)
        
        
    end
    waitbar(iteration/(iteration_last-1),f,'Loading data for figures','Name','Reloading previously obtained data');
end
maxS=max(mag_spec);
limS_updated=limS;
while limS_updated<=maxS
    limS_updated=limS_updated*limSstretcher;
end
if figure_analysis==3
    E=0;ac_t=[0 0 0];
end

delete(f)

end

% function ac_t_current=min_max_ave_allen_cahn(allen_cahn)
% ac_t_current=[max(max(allen_cahn)) min(min(allen_cahn)) mean(mean(allen_cahn))];
% end