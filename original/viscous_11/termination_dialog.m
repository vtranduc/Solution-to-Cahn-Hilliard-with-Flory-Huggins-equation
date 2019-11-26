function termination_dialog()
button = questdlg('Simulation has been cancelled.','Simulation is shutting down.','Ok','Ok');
switch button
    case 'Ok'
        return
end
end