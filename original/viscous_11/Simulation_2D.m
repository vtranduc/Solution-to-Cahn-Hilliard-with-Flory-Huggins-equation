function Simulation_2D

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

diff=5000;

T=0.7;

% ci = 0.62; %Overall initial concentrations
ci=0.7
% ci=[0.5,0.6];
% ci=0.7
include_fluc=1;
ci_fluc = 0.01; %Initial fluctuation

nr_tol=1.0e-6;
nr_max_iteration=20;

time_out=3600*(1/6);

entropy=1;
T_theta=1;

fps=10;

dt_down=0.5;
dt_up=1.2;
dt_min=1.0e-20;
dt_ideal=5.0e-8;
bypass_tol=5;

max_frame=20;

xlen=1;ylen=1;

alpha=200;
beta=8;

frame_size=[1920 1080];
frame_pos=[0 0];

limS=1; %Guessed value. This will be adjusted as necessary
limShigher=1;
limSlower=0.01;
limSstretcher=10;
limScompressor=0.1;

fc=15; %Put fc=0 if you wanna use Nyquist frequency

num_str_fac_1d=2;

limS1d=1;
limS1dhigher=1;
limS1dlower=0.01;
limSstretcher1d=10;
limScompressor1d=0.1;

nT_phase=1000000000;
T_min_phase=0.001;

n1=1;
n2=10;

ifig_phase=9;

proceed_request=0; % Yes=1

%Phase diagram specifications
nT=100;
tol=1.0e-6;
parallel_computing=0;
nInitialGuessPts=100;
guess_accuracy=10^-3;

illustration=1; %1=yes, 0=no

% === Validate user inputs ================================================

%COMPLETE IT LATER IN ANOTHER FUNCTION!

%CHECK TO MAKE SURE CI+/-CI_FLUC IS WITHIN 0 AND 1!!! 

%Entropy, n1, and T_theta must be 1!!!

%T_min should be larger than zero. Otherwise, the end value (T=0) will yield corresponding spinodal point as NaN. Should be as close to 0 as possible

% === Illustrate phase diagram ============================================

fprintf('Generating phase diagram...\n')

if illustration==1
    [xspinodal,yspinodal,xbinodal,ybinodal]=spinodal_binodal(...
        n1,n2,entropy,T_theta,nT,tol,parallel_computing,nInitialGuessPts,...
        guess_accuracy);
    figure(ifig_phase)
    illustrate_simulation(xspinodal,yspinodal,xbinodal,ybinodal,T,ci)
end

if proceed_request==1
    proceed=input('Proceed? Yes=1+ENTER or Yes=ENTER, No=Any other numerical+ENTER or string+ENTER\n');
    if proceed~=1
        disp('The simulation has been cancelled by the user')
        return
    else
        disp('The simulation is to proceed')
    end
end

% === Specify constants which will be used throughout =====================

w=[5/18 4/9 5/18];
gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

nx=nex+1;ny=ney+1;
ne=nex*ney;
n=nx*ny;
nfour=4*n;

%===========TESTING================================

% nnz_=nnzSJ(nx,ny);
% % 
% % A=ones(1,nnz_);
% % b=ones(1,nnz_);
% % c=ones(1,nnz_);
% % d=ones(1,nnz_);
% % e=ones(1,nnz_);
% % f=ones(1,nnz_);
% % g=ones(1,nnz_);
% % h=ones(1,nnz_);
% % i=ones(1,nnz_);
% % j=ones(1,nnz_);
% % k=ones(1,nnz_);
% 
% lll=ones(20,nnz_);
% kk=ones(100,nfour);
% 
% [edge_xdir,edge_ydir]=edgy(10,ny,n);
% 
% [irow,icol]=setUpSparse(nnz_,ny,n);
% 
% [irow' icol'];
% 
% domain=linspace(0,1,1000);
% 
% plot(domain,basis(domain,1,1,0))
% grid on
% 
% weights=generateWeights(); %This has not been verified yet
% 
% display('Testing up to here')
% return

%======================================================

ll=2;li=4;lu=2+4*ney;
bl=3;bi=4*ny;bu=3+4*nex*ny;
tl=3+4*ney;ti=4*ny;tu=3+4*ney+4*nex*ny;
rl=2+4*nex*ny;ri=4;ru=2+4*ney+4*nex*ny;

if fc==0
    fc=1/(2*min([xlen/nex ylen/ney]));
end

frequency_domain=linspace(0,fc,alpha);

if num_str_fac_1d>nx
    error('num_str_fac_1d must be smaller than or equal to nx!')
end
isf1d=zeros(1,num_str_fac_1d);
spacer=nx/(num_str_fac_1d-1);
curr_pos=0;
for i=2:1:num_str_fac_1d-1
    curr_pos=curr_pos+spacer;
    isf1d(i)=round(curr_pos);
end
isf1d(1)=1;
isf1d(num_str_fac_1d)=nx;

% === initialize time =====================================================

time=0.0;

% === Evaluate value of chi ===============================================

chi=0.5-entropy*(1-T_theta/T)

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

% co=co_linear_gradient_x(low_c,high_c,ci_fluc,nx,ny);
co=generate_co_2D(ci,nx,ny,include_fluc,ci_fluc);
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
freq_domain=linspace(0,fc,alpha);
structure_factor_1d=structure_factor_1D_multiple(c_nodal,...
    freq_domain,isf1d);


% === Store & Adjust the first structure factor data ======================

maxS=max(structure_factor);
while limS<=maxS
    limS=limS*limSstretcher;
end
logS=[time log(maxS)];

maxS1d_mult=max(structure_factor_1d,[],2);
maxS1d=max(maxS1d_mult);
while limS1d<=maxS1d
    limS1d=limS1d*limSstretcher1d;
end
logS1d=log(maxS1d_mult);

% === Set up for plotting =================================================

fri=0;
if illustration==1
    fig=figure(1);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    minz=0.0;maxz=1.0;
    

% === Plot the concentration in second step, which is qual to first one ===

    fri=fri+1;
    visual(fri)=visual2D_gradient(c,nx,ny,x_coord,y_coord,minz,maxz,time,...
        frequency_domain,structure_factor,ci_ave,limS,logS,...
        structure_factor_1d,freq_domain,logS1d,limS1d,...
        xspinodal,yspinodal,xbinodal,ybinodal,T,ci);
end

% === Prepare to enter main loop ==========================================

enopm=zeros(1,4);
enop=zeros(1,16);
cp=zeros(1,nfour);
bypass=0;

% === Clear variables that will no longer be needed

%CLEAR ALL VARIABLES THAT WILL NO LONGER BE NEEDED!!!

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
        
%         %-----Testing-----------
% toc
%         tester111=compute_sf(...
%             nfour,weight_adjuster(generateWeights(),dx,dy),ny,n,chi*ones(1,nx),c,co,diff,T*ones(1,nx),n1,n2,dt,dx,dy);
% %         nfour,weights,ny,n,chi,c,co,diff,T,n1,n2,dt,dx,dy
%         myDiff=abs(sf+tester111);
%         for ijk=1:1:nfour
%             if abs(myDiff(ijk))>=0.0001
%                 myDiff(ijk)
%                 error('bad diff')
%             end
%         end
%         [sf' tester111' myDiff']
%         max(myDiff)
%         min(myDiff)
%         toc
%         warning('dfafdsa')
%         nnz_=nnzSJ(nx,ny);
%         [irow,icol]=setUpSparse(nnz_,ny,n);
%         weights=generateWeights();
%         weights=weight_adjuster(weights,dx,dy);
%         [adjusted_chi,adjusted_diffT]=chi_diffT_adjuster(T*ones(1,nx),dx,entropy,T_theta,diff);
%         tester222=compute_sj(nnz_,irow,icol,ny,n,c,weights,ney,dt,...
%             adjusted_chi,n1,n2,adjusted_diffT,dx*dy);
%         tester333=zeros(nfour,nfour);
%         for iiiii=1:1:length(irow)
%             tester333(irow(iiiii),icol(iiiii))=tester222(iiiii);
%         end
%         sj;
%         reshapesj=reshape(sj,1,nfour*nfour);
%         reshapetester333=reshape(tester333,1,nfour^2);
%         tester555=abs(reshapesj-reshapetester333);
%         [reshapesj' reshapetester333' tester555']
%         max(tester555)
%         min(tester555)
%         warning('Systematic love')
%         %-----------------------

        [sf,sj]=bcsfsj(sf,sj,nx,ny,...
            0,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);
        
        
        % --- Testing Boundary conditions ---------
%         disp('Testing boundary conditions')
%         tester1000=tester222;
%         tester222=test23(tester222,nx,ny,nnz_,irow,icol,sj);
%         tester333=zeros(nfour,nfour);
%         tester1001=zeros(nfour,nfour);
%         for iiiii=1:1:length(irow)
%             tester333(irow(iiiii),icol(iiiii))=tester222(iiiii);
%             tester1001(irow(iiiii),icol(iiiii))=tester1000(iiiii);
%         end
%         
%         tester2000=abs(tester333-sj);
%         
%         for iirow=1:1:nfour
%             for iicol=1:1:nfour
%                 if tester2000(iirow,iicol)>1
%                     
%                     iirow
%                     iicol
%                     sj(iirow,iicol)
%                     tester2000(iirow,iicol)
%                     warning('let us try to find it in sparse array')
%                     for iaadfadf=1:1:length(irow)
%                         if irow(iaadfadf)==iirow && icol(iaadfadf)==iicol
%                             disp('Found!')
%                             break
%                         end
%                     end
%                     iaadfadf
%                     tester222(iaadfadf)
%                     sj(iirow,iicol)
%                     irow(1);
%                     icol(1);
%                     error('too big lala')
%                 end
%             end
%         end
%         
%         reshapesj=reshape(sj,1,nfour*nfour);
%         reshapetester333=reshape(tester333,1,nfour^2);
%         tester555=abs(reshapesj-reshapetester333);
%         tester1001=reshape(tester1001,1,nfour^2);
%         
%         [reshapesj' reshapetester333' tester1001' tester555']
%         max(tester555)
%         min(tester555)
%         disp('dangerous')
%         max(max(tester2000))
        
        
%         disp('New method')
% %         test=nnzsj(500,500)
%         [iZeros,iOnes]=get_bc_keyVals(nx,ny,nnz_,irow,icol);
%         
%         [tester1212,tester333]=parallel_bc(tester222,iZeros,iOnes,...
%             tester111,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu);
%         
%         tester555=zeros(nfour,nfour);
%         for owaranai=1:1:length(irow)
%             tester555(irow(owaranai),icol(owaranai))=tester333(owaranai);
%         end
%         tester666=reshape(tester555,1,nfour^2);
%         tester777=reshape(sj,1,nfour^2);
%         
%         tester888=abs(tester777-tester666);
%         
%         [max(tester888) min(tester888)]
%         
%         tester1212=-tester1212
%         tester999=abs(tester1212-sf);
%         
%         
%         [sf' tester999' tester999']
%         
%         [max(tester999) min(tester999)]
%         
%         error('Systematic love')
        %------------------------------------------
        
% det(sj)
%         sf'
%         error('dfdasdf just stop')

% === Create compact matrix ===============================================
        fprintf('Compressing')
        toc
%         [rowi,coli,val]=compress_sj(sj,nfour,nfour);
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
    
    structure_factor_1d=structure_factor_1D_multiple(c_nodal,...
        freq_domain,isf1d);
    
% === Analyze structure factor ============================================
    
    fri=fri+1;
    
    if illustration==1
        maxS=max(structure_factor);
        while limS<=maxS
            limS=limS*limSstretcher;
        end
        logS(fri,1:1:2)=[time log(maxS)];

        maxS1d_mult=max(structure_factor_1d,[],2);
        maxS1d=max(maxS1d_mult);
        while limS1d<=maxS1d
            limS1d=limS1d*limSstretcher1d;
        end
        logS1d(1:1:num_str_fac_1d,fri)=log(maxS1d_mult);
%     end
    
    
% === Visualize the result ================================================

%     if illustration==1
        visual(fri)=visual2D_gradient(c,nx,ny,x_coord,y_coord,minz,maxz,time,...
            frequency_domain,structure_factor,ci_ave,limS,logS,...
            structure_factor_1d,freq_domain,logS1d,limS1d,...
            xspinodal,yspinodal,xbinodal,ybinodal,T,ci);
    end
    
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
    
    % WRONG!
    if time>obs_t
        time=obs_t;
        dt=obs_t-(time-dt);
    end
end

% === Export the movie ====================================================
toc
if illustration==1
    elapsed_=cputime-counter_init;
    if length(ci)==1
        vid_title=sprintf('2D nonlinear Cahn Hilliard,%dx%d,ci=%1.3f,fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci_ave,ci_fluc*include_fluc,diff,T,elapsed_);
    elseif length(ci)==2
        vid_title=sprintf('2D nonlinear Cahn Hilliard,%dx%d,ci=[%1.2f,%1.2f],fluc=%1.4f,diff=%d,T=%0.3f,elapsed time=%1.3fs.avi',nex,ney,ci(1),ci(2),ci_fluc*include_fluc,diff,T,elapsed_);
    end
    video=VideoWriter(vid_title, 'Uncompressed AVI');
    video.FrameRate=fps;
    open(video)
    writeVideo(video,visual);
    close(video)


% === Identify transition time ============================================

    figure(2)
    plot(logS(:,1),logS(:,2),'x')
    grid on

% === Conclusion ==========================================================
% dlmwrite('fourier.txt',c)


    fig=figure(3);
    set(fig,'Position',[frame_pos(1) frame_pos(2) frame_size(1) frame_size(2)])
    movie(gcf,visual,1)
end
toc
end