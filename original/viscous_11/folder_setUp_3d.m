function output_path=folder_setUp_3d(export_data)

output_path=1;
while 1
    if ~exist([pwd '/Results_3d/Run_' num2str(output_path) '/'],'dir')
        break
    end
    output_path=output_path+1;
end
output_path=[pwd '/Results_3d/Run_' num2str(output_path) '/'];
mkdir(output_path);
if export_data==1
    mkdir([output_path 'specification/'])
    mkdir([output_path 'time/'])
    mkdir([output_path 'concentration/'])
end

end