function ac_t=test82(ac_max,ac_min,ac_ave,t)
% ac_max=1;
% ac_min=0;
% ac_ave=0.5;
% t=0;
if t==0 && exist('nozomi.dat','file')
    delete('nozomi.dat')
end
if ~exist('nozomi.dat','file')
    ac_t=[ac_max,ac_min,ac_ave];
    dlmwrite('nozomi.dat',ac_t);
    return
end
ac_t=dlmread('nozomi.dat');
[m,~]=size(ac_t);
ac_t(m+1,:)=[ac_max,ac_min,ac_ave];
dlmwrite('nozomi.dat',ac_t);
end