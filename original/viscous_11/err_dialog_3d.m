function err_dialog_3d(message)
button = questdlg(message,'Error','Ok','Ok');
switch button
    case 'Ok'
        return
end
end