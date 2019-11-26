function test65

clear
clc



eX=[0 2 0 2 0 2 0 2];
eY=[-2 -1 2 1 -2 -1 2 1];
eZ=[-1 -1 -1 -1 1 1 1 1];

alpha=0;beta=0.2;gamma=0.1;


terms_to_be_obtained=[1 2 3];

take=get_fg_dervs_3d(alpha,beta,gamma,[0 0 0],[0 1 0],eX,eY,eZ,terms_to_be_obtained);

n=[-2 -0.6 0];

take'

dot(n,take)

end