function nonlinear_Cahn_Hilliard

clear
clc

% ======================= Start the counter ===============================

tic
counter_init=cputime;

% ===================== User defined variables ============================

nex=20;
ney=20;

dto=1.0e-5;
dt=1.0e-5;

obs_t=2.0e-1; %Observation time

diff=1800;

T=0.6;

ci = 0.7; %Overall initial concentrations
ci_fluc = 0.01; %Initial fluctuation

nr_tol=1.0e-6;
nr_max_iteration=20;

time_out=3600*0.5;

entropy=1;
T_theta=1;

fps=10;

dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-20;
dt_ideal=5.0e-8;
bypass_tol=5;

max_frame=100;

xlen=1;ylen=1;

alpha=200;
beta=8;

frame_size=[1920 1080];
frame_pos=[0 0];

limS=16; %Guessed value. This is adjusted as program progresses

fc=12; %Put fc=0 if you wanna use Nyquist frequency

% === Specify constants which will be used throughout =====================

w=[0.27778 0.4444 0.27778];
gp=[0.1127 0.5 0.8873];

nx=nex+1;ny=ney+1;
ne=nex*ney;
n=nx*ny;
nfour=4*n;

ll=2;li=4;lu=2+4*ney;
bl=3;bi=4*ny;bu=3+4*nex*ny;
tl=3+4*ney;ti=4*ny;tu=3+4*ney+4*nex*ny;
rl=2+4*nex*ny;ri=4;ru=2+4*ney+4*nex*ny;

if fc==0
    fc=1/(2*min([xlen/nex ylen/ney]));
end

frequency_domain=linspace(0,fc,alpha);

% === initialize time =====================================================

time=0.0;

% === Evaluate value of chi ===============================================

chi=0.5-entropy*(1-T_theta/T);

% === Determine the coordinate at each node ===============================

x=zeros(1,n); y=zeros(1,n);

for i=1:n
    x(i)=(xlen/nex)*floor((i-1)/ny);
end
for i=1:ny
    for j=1:nx
        y(i+(j-1)*ny)=(ylen/ney)*(i-1);
    end
end

x_coord=linspace(0,1,nx);
y_coord=linspace(0,1,ny);

% === Generate matrices that help to characterize each element uniquely ===

nop=zeros(ne,16);
for i=1:ne
    nop(i,1)=4*floor((i-1)/ney)+4*(i-1)+1;
    nop(i,9)=nop(i,1)+4*ny;
    for j=1:7
        nop(i,1+j)=nop(i,1)+j;
        nop(i,9+j)=nop(i,9)+j;
    end
end
nopm=zeros(i,4);
for i=1:ne
    nopm(i,1)=floor((i-1)/ney)+i;
    nopm(i,2)=nopm(i,1)+1;
    nopm(i,3)=nopm(i,1)+ny;
    nopm(i,4)=nopm(i,3)+1;
end

% === Randomize and apply initial conditions ==============================

co = zeros(1,nfour);

for i=1:1:n
    co((i-1)*4+1)=ci+(ci_fluc*(2*rand(1,1)-1));
end

c_nodal=extract_nodal_weights_2D(co,nx,ny);
ci_ave=mean(mean(c_nodal));

ci_diff=ci-ci_ave;

for i=1:1:n
    j=(i-1)*4+1;
    co(j)=co(j)+ci_diff;
end

co=icc(co,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);

% === Set up Residual vector and Jacobian matrix ==========================


fprintf('Setting up')
toc
% special_zero=1.0e-100;

% === Take the very first time step =======================================

%time=time+dto;
c=co;

% === Analyze initial concentration =======================================

c_nodal=extract_nodal_weights_2D(c,nx,ny);
ci_ave=mean(mean(c_nodal));

% === Compute first structural factor =====================================

structure_factor=structure_factor_2D(alpha,beta,fc,c_nodal,ci_ave);

% === Store & Adjust the first structure factor data ======================

maxS=max(structure_factor);
while limS<=maxS
    limS=limS*2;
end
logS=[time log(maxS)];

% === Set up for plotting =================================================

fig=figure(1);
set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
minz=0.0;maxz=1.0;
fri=0;

% === Plot the concentration in second step, which is qual to first one ===

fri=fri+1;
visual(fri)=visual2D(c_nodal,x_coord,y_coord,minz,maxz,time,...
    frequency_domain,structure_factor,ci_ave,limS,...
    logS);

% === Prepare to enter main loop ==========================================

enopm=zeros(1,4);
enop=zeros(1,16);
cp=zeros(1,nfour);
bypass=0;

% === Enter timer loop ====================================================
fprintf('Enter main cycle')
%for rahmat=1:1:20
while 1

% === Update the concentration data =======================================
    
    for i=1:1:nfour
        cp(i)=c(i)+dt*((c(i)-co(i))/dto);
    end
    coo=co;
    co=c;
    c=cp;
   
% === Apply BC to newly predicted concentration ===========================

    c=bc(c,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);

% === Enter Newton-Raphson iterations =====================================

    jk=0;
    err = nr_tol+1;
    
	while err > nr_tol

% === Reset residual vector and Jacobian matrix ===========================

        sf=zeros(1,nfour);
        sj=zeros(nfour,nfour);

% === Start looping through each element ==================================
        
        fprintf('Filling in sf sj')
        for e=1:1:ne
            
% === Characterize the current element ====================================

            for i=1:1:4
                enopm(i)=nopm(e,i);
            end
            for i=1:1:16
                enop(i)=nop(e,i);
            end
            dx=x(enopm(3))-x(enopm(1));
            dy=y(enopm(2))-y(enopm(1));
% [dx dy]
% === Shift the weights of two basis functions ============================
            
            for w1=1:1:3
                for w2=1:1:3

% === Characterize the current element's basis accordingly ================

                    [phi,phix,phiy,phixx,phiyy,~]=...
                        tfunct(gp(w1),gp(w2),dx,dy);

% === Determine previous and current absolute values, and associated slopes
% at current particular pint ==============================================

                    con=0.0;
                    cono=0.0;
                    conx=0.0;
                    cony=0.0;
                    conxx=0.0;
                    conyy=0.0;

                    for i=1:1:16
                        con = con + c(enop(i)) * phi(i);
                        cono = cono + co(enop(i)) * phi(i);
                        conx = conx + c(enop(i)) * phix(i);
                        cony = cony + c(enop(i)) * phiy(i);
                        conxx = conxx + c(enop(i)) * phixx(i);
                        conyy = conyy + c(enop(i)) * phiyy(i);
                    end
                    
% === Change in concentration over time step is estimated =================

                    cont=(con-cono)/dt;
                    
% === Fill in residual vector =============================================

                    for i=1:1:16
                        sf(enop(i))=sf(enop(i))-w(w1)*w(w2)*dx*dy*...
                            (cont*phi(i)-diff*T*phi(i)*(-1.0/(con^2.0)...
                            +(1./10.)*(1-con)^-2.0)*(conx^2.0+...
                            cony^2.0)-diff*T*phi(i)*...
                            (1/con+(1./10.)*(1./(1.-con))-2*chi)*...
                            (conxx+conyy)+(conxx+conyy)*...
                            (phixx(i)+phiyy(i)));
                        
% === Fill in the Jacobian matrix =========================================

                        for j=1:1:16
                            sj(enop(i),enop(j))=sj(enop(i),enop(j))+w(w1)...
                                *w(w2)*dx*dy*...
                                (phi(i)*phi(j)/dt-diff*T*phi(i)*((2*phi(j)...
                                /(con^3.0)+2*phi(j)/(10*(1-con)^3))*...
                                (conx^2.0+cony^2.0)+(-1/con^2.0+...
                                1/(10*(1-con)^2.0))...
                                *(2*conx*phi(j)+2*cony*phi(j)))-diff*T*phi(i)...
                                *((-phi(j)/con^2.0+phi(j)/(10*(1-con)^2))*...
                                (conxx+conyy)+(1/con+1/(10*(1-con))-2*chi)...
                                *(phixx(j)+phiyy(j)))+(phixx(j)+phiyy(j))...
                                *(phixx(i)+phiyy(i)));
                            
                        end
                    end
                end
            end
        end
% det(sj)
% === Apply boundary conditions to sf and sj ==============================
        fprintf('Apply bc')
        [sf,sj]=bcsfsj(sf,sj,nx,ny,...
            0,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);
% det(sj)

% === Create compact matrix ===============================================
        fprintf('Compressing')
        toc
        %[rowi,coli,val]=compress_sj(sj,nfour,nfour);
        [rowi,coli,val]=find(sj);
        toc
        %[rowi',coli',val']
        fprintf('Sparsing')
        sj=sparse(rowi,coli,val);
        toc
        
% === Carry matrix Division ===============================================
        fprintf('Carry out division')
        c_= sj \ sf';
%         error('just stop')
% === Update the solution =================================================

        c = c + c_';
        
% === Evaluate the error ==================================================
        fprintf('error eval')
        err = sqrt(sum(c_.^2.0))
        
% === Assert reasonable convergence =======================================

        jk = jk +1;
        if jk >=nr_max_iteration || (jk>1 && err>=erro)
            %error('Newton-Raphson is not converging')
            fprintf('Newton-Raphson is not converging')
            jk=-1;
            break
        else
            erro=err;
        end
	end
    
% === Resolve divergence ==================================================
    
    if jk==-1
        bypass=bypass+1;
        c=co;
        co=coo;
        if bypass>bypass_tol
            fprintf('Too much bypassing!')
            break
        else
            dt=dt*dt_down;
            if dt<dt_min
                fprintf('Time step is too small!')
                break
            else
                continue
            end
        end
    elseif bypass~=0
        bypass=0;
    end
    
% === Increment time step =================================================

    time=time+dt % Time for last loop, not next!
    
% === Extract nodal weights ===============================================

    fprintf('Extracting nodal weights')
    c_nodal=extract_nodal_weights_2D(c,nx,ny);
    
% === Compute structure factor ============================================

    fprintf('Computing structure factor')
    structure_factor=structure_factor_2D(alpha,beta,fc,c_nodal,ci_ave);
    
% === Analyze structure factor ============================================
    
    fri=fri+1;
    maxS=max(structure_factor);
    while limS<=maxS
        limS=limS*2;
    end
    logS(fri,1:1:2)=[time log(maxS)];
    
% === Visualize the result ================================================

    visual(fri)=visual2D(c_nodal,x_coord,y_coord,minz,maxz,time,...
        frequency_domain,structure_factor,ci_ave,limS,...
        logS);
    if fri==max_frame
        fprintf('Max frame reached!')
        break
    end
    
% === Store structure factor data =========================================

    logS(fri,1:1:2)=[time log(max(structure_factor))];
    
% === Time out check point ================================================

    if time_out<cputime-counter_init
        fprintf('Time out!')
        break
    end
    
% === Prepare for next time loop ==========================================

    if time>=obs_t
        fprintf('Observation time goal reached!')
        break
    end
    
    dtoo=dto;
    dto=dt;
    
    if jk>=10
        dt=dt*dt_down;
        if dt<dt_min
            fprintf(' Time step is too small! ')
            break
        end
    elseif jk<=3 || dt<dt_ideal
        dt=dt*dt_up;
    end
    
    if time>obs_t
        time=obs_t;
        dt=obs_t-(time-dt);
    end
end

% === Export the movie ====================================================

elapsed_=cputime-counter_init;
vid_title=sprintf('2D nonlinear Cahn Hilliard, %d x %d, ci=%1.5f, diff=%d,T=%d, elapsed time=%d .avi',nex,ney,ci_ave,diff,T,elapsed_);
video=VideoWriter(vid_title, 'Uncompressed AVI');
video.FrameRate=fps;
open(video)
writeVideo(video,visual);
close(video)

% === Identify transition time ============================================

figure(2)
plot(logS(:,1),logS(:,2),'x')

% === Conclusion ==========================================================
% dlmwrite('fourier.txt',c)
toc

fig=figure(3);
set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
movie(gcf,visual,1)

end