function test25

clear
clc

A=ones(2,2)

currentFolder = pwd

folder_name='testOutput'

trial=2;

output_path=[pwd '/testOutput_' num2str(trial) '/']

ifig=5;

conc_path=[output_path 'conc_' num2str(ifig) '.txt']

mkdir('testOutput_2')


dlmwrite(conc_path,A)

if exist('testOutput_2','dir')
    disp('heel')
end

end