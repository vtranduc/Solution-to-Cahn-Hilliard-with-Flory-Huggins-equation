function proceed=proceed_dialog()
button = questdlg('Ready to proceed?','Phase diagram evaluation','Yes','No','No');
switch button
    case 'Yes'
        disp('Simulation is to proceed')
        proceed=1;
        return
    case 'No'
        disp('Simulation has been cancelled by the user')
        proceed=0;
end

end