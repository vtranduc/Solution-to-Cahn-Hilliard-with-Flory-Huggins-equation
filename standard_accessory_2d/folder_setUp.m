function output_path=folder_setUp(export_data,export_figure,export_spec,k_specific)

output_path=1;
while 1
    if ~exist([pwd '/Results/Run_' num2str(output_path) '/'],'dir')
        break
    end
    output_path=output_path+1;
end
output_path=[pwd '/Results/Run_' num2str(output_path) '/'];
mkdir(output_path);

if export_data==1
    mkdir([output_path 'specification/'])
    mkdir([output_path 'time/'])
    mkdir([output_path 'concentration/'])
end
if export_figure==1
    mkdir([output_path 'figures/'])
    if export_spec==1
        mkdir([output_path 'figures/simulation_range/'])
        mkdir([output_path 'figures/profile/'])
        mkdir([output_path 'figures/phases/'])
        mkdir([output_path 'figures/fourier_transform/'])
        mkdir([output_path 'figures/stage_view/'])
        mkdir([output_path 'figures/energy/'])
        mkdir([output_path 'figures/max_c/'])
        mkdir([output_path 'figures/min_c/'])
    elseif export_spec==2
        mkdir([output_path 'figures/potential/'])
    elseif export_spec==3
        mkdir([output_path 'figures/k=' num2str(k_specific) '/'])
    end
end

end