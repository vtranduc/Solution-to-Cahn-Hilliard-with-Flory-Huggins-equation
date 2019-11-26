function second_trial

clear
clc

tic

%Page 1

W=[0.27778 0.4444 0.27778];
GP=[0.1127 0.5 0.8873];

NEX = 10;
NEY = 10;
NX = NEX + 1;
NY = NEY + 1;
NE = NEX * NEY;
N = NX * NY;
NFOUR = 4 * N;
TIME = 0.0D0;
format short
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
%-----Good for 100x100 elememts----------
%DTO = 5.0D-10;
%DT = 5.0D-10;
%max_time=1.0D-5;
%-----Perfect for 10x10 elememts----------
%{
DTO = 1.0D-5;
DT = 1.0D-5;
max_time=2.0D-3;
max_time_step=1.0D-4;
%}
%-----------------------------------------

%--------Good for 10x10 nonlinear Cahn-Hilliard--------
%{
DTO = 1.0D-5;
DT = 1.0D-5;
max_time=2.0D-4;
%}
%------------------------------------------------

%----VERIFIED ~90 MIN RUNTIME 100x100 nonlinear Cahn-Hilliard--------
%{
DTO = 1.0D-5;
DT = 1.0D-5;
max_time=2.0D-4;
%}
%------------------------------------------------

%------------------------------------------------
%%{
DTO = 1.0D-5;
DT = 1.0D-5;
max_time=2.0D-4;

DT_catchup=5.0D-6;
catchup_pace=1.5;
%%}
%------------------------------------------------

%-----Almost perfect for 100x100 elememts----------
%{
DTO = 1.0D-11;
DT = 1.0D-11;
max_time=2.5D-3;
max_time_step=1.0D-4;
%}
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
DIFF = 6.0D2;
%=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

frame_index=0;
bypass_counter=0;
%=======================================================================2

X=zeros(1,N);
Y=zeros(1,N);

for I=1:N
    X(I)=(1/NEX)*floor((I-1)/NY);
end
for I=1:NY
    for J=1:NX
        Y(I+(J-1)*NY)=(1/NEY)*(I-1);
    end
end
NOP=zeros(NE,16);
for I=1:NE
    NOP(I,1)=4*floor((I-1)/NEY)+4*(I-1)+1;
    NOP(I,9)=NOP(I,1)+4*NY;
    for J=1:7
        NOP(I,1+J)=NOP(I,1)+J;
        NOP(I,9+J)=NOP(I,9)+J;
    end
end
NOPM=zeros(I,4);
for I=1:NE
    NOPM(I,1)=floor((I-1)/NEY)+I;
    NOPM(I,2)=NOPM(I,1)+1;
    NOPM(I,3)=NOPM(I,1)+NY;
    NOPM(I,4)=NOPM(I,3)+1;
end
CO=zeros(1,NFOUR);
%RANDOMIZE IC==============================================================
for I=1:NFOUR
    if rand(1)<=0.5
        mult=-1;
    else
        mult=1;
    end
    CO(I)=0.7+mult*rand(1)*0.01;
end

%==========================================================================
K=0;
for J=1:4:1+4*NEY
    K=K+1;
    K1=K;
    for I=J:4*NY:J+4*NEX*NY
        K1=K1+NY;
    end
end
TIME=TIME+DTO;

%=======================================================================3
C=zeros(1,NFOUR);
for J = 1 : NFOUR
    C(J) = CO(J);
end
%Left
for J=2:4:2+4*NEY
    J1 = J + 1 ;
    J2 = J + 2;
    C(J) = 0.000;
    C(J1) = 0.000;
    C(J2) = 0.000;
end
%Bottom
for J=2:4*NY:2+4*NEX*NY
    J1=J+1;
    J2 = J + 2;
    C(J) = 0.000;
    C(J1) = 0.000;
    C(J2) = 0.000;
end
%Top
for J = 2 + 4 * NEY:4 * NY: 2 + 4 * NEY + 4 * NEX * NY 
    J1 = J + 1;
    J2 = J + 2;
    C(J) = 0.000;
    C(J1) = 0.000;
    C(J2) = 0.000;
end
%Right
for J = 2 + 4 * NEX * NY: 4: 2 + 4 * NEY + 4 * NEX * NY
    J1 = J + 1;
    J2 = J + 2;
    C(J) = 0.000;
    C(J1) = 0.000;
    C(J2) = 0.000;
end

%{
for J = 1 : NFOUR
    C(J) = CO(J);
end
%Left
for J=2:4:2+4*NEY
    %J1 = J + 1 ;
    %J2 = J + 2;
    C(J) = 0.000;
    %C(J1) = 0.000;
    %C(J2) = 0.000;
end
%Bottom
for J=2:4*NY:2+4*NEX*NY
    J1=J+1;
    %J2 = J + 2;
    %C(J) = 0.000;
    C(J1) = 0.000;
    %C(J2) = 0.000;
end
%Top
for J = 2 + 4 * NEY:4 * NY: 2 + 4 * NEY + 4 * NEX * NY 
    J1 = J + 1;
    %J2 = J + 2;
    %C(J) = 0.000;
    C(J1) = 0.000;
    %C(J2) = 0.000;
end
%Right
for J = 2 + 4 * NEX * NY: 4: 2 + 4 * NEY + 4 * NEX * NY
    %J1 = J + 1;
    %J2 = J + 2;
    C(J) = 0.000;
    %C(J1) = 0.000;
    %C(J2) = 0.000;
end
%}

%PLOT INITIAL===========================
x_coord=linspace(0,1,NX);
y_coord=linspace(0,1,NY);
%===================
solution=zeros(1,NX*NY);
for i=1:NX*NY
    solution(i)=C((i-1)*4+1);
end

%max_z=max(solution)*1.;
%min_z=min(solution)*1.;

max_z=1;
min_z=0;

max_bar=0.51;
min_bar=0.49;

%myave=mean(mean(mean(CO)));
%max_bar=myave+0.001;
%min_bar=myave-0.001;

%max_bar=0.501;
%min_bar=0.499;

solution=reshape(solution,[NY,NX]);
frame_index=frame_index+1;

hFig = figure(1);
pos_size=[0 200 1280 720];
set(hFig, 'Position', pos_size)

%INITIAL PLOT========================================================
subplot(1,2,1)
surf(x_coord,y_coord,solution);
colorbar('southoutside')
%caxis([min_bar max_bar]);
axis([0 1 0 1 min_z max_z])
xlabel('x'),ylabel('y'),zlabel('z')
subplot(1,2,2)
contourf(x_coord,y_coord,solution);
%set(hC,'LineStyle','none');
colorbar
%caxis([min_bar max_bar]);
suptitle('Solution to 2D nonlinear Cahn-Hilliard in unstable region')
myClip(frame_index)=getframe(gcf);
%======================================================================

CP=zeros(1,NFOUR);
COO=zeros(1,NFOUR);

bypassed=0;

CHI=1.1666666666667;

toc

display('Start main loop')

while 1
%for i=1:1:1
    % PREDICT VALUES FOR NEXT TIME STEP
    
    display('Prediction starts')
    
    for J = 1: NFOUR
        CP(J) = C(J) + DT * ( (C(J)-CO(J) ) / DTO);
        COO(J) = CO(J);
        CO(J) = C(J);
        C(J) = CP(J);
    end
    for J = 2: 4:2 + 4 * NEY
        J1 = J + 1;
        J2 = J + 2;
        C(J) = 0.000;
        C(J1) = 0.000;
        C(J2) = 0.000;
    end
    for J = 2:4 * NY: 2 + 4 * NEX * NY
        J1 = J + 1;
        J2 = J + 2;
        C(J) = 0.000;
        C(J1) = 0.000;
        C(J2) = 0.000;
    end

    %=======================================================================4

    for J=2+4*NEY:4*NY:2+4*NEY+4*NEX*NY
        J1=J+1 ;
        J2=J+2;
        C(J) = 0.000;
        C(J1) = 0.000;
        C(J2) = 0.000;
    end
    for J=2+4*NEX*NY:4:2+4*NEY+4*NEX*NY
        J1=J+1;
        J2=J+2;
        C(J) = 0.000;
        C(J1) = 0.000;
        C(J2) = 0.000;
    end

    INOP=zeros(1,16);
    INOPM=zeros(1,4);
    
    display('Newton-Raphson starts')
    
    JK = 0;

    %while 1
    for i=1:2
        try
            SJ=zeros(NFOUR,NFOUR);
            SF=zeros(1,NFOUR);
        catch
            SJ=sparse(NFOUR,NFOUR);
            SF=sparse(1,NFOUR);
        end
        % ENTER THE 100 D0 LOOP
        display('Start filling in vectors')
        for I = 1: NE
            for J = 1: 4
                INOPM(J) = NOPM(I, J);
            end
            for J = 1: 16
                INOP(J) = NOP(I, J);
            end
            DX = X(INOPM(3)) - X(INOPM(1));
            DY = Y(INOPM(2)) - Y(INOPM(1));
            for J = 1: 3
                for K = 1: 3
                    [PHI,PHIX,PHIY,PHIXX,PHIYY,PHIXY]=TFUNCT(GP(J),GP(K),DX,DY);
                    CON = 0.000;
                    CONO = 0.000;
                    CONX = 0.000;
                    CONY = 0.000;
                    CONXX = 0.000;
                    CONYY = 0.000;
                    for INTR = 1: 16
                        CON = CON + C(INOP(INTR)) * PHI(INTR);
                        CONO = CONO + CO(INOP(INTR)) * PHI(INTR);
                        CONX = CONX + C(INOP(INTR)) * PHIX(INTR);
                        CONY = CONY + C(INOP(INTR)) * PHIY(INTR);
                        CONXX = CONXX + C(INOP(INTR)) * PHIXX(INTR);
                        CONYY = CONYY + C(INOP(INTR)) * PHIYY(INTR);
                    end
                    DCON = (CON - CONO) / DT;
                    %=======================================================================5
                    for L = 1: 16
                        %{
                        SF(INOP(L))=SF(INOP(L))-W(J)*W(K)*DX*DY*...
                            (DCON*PHI(L)-DIFF*PHI(L)*(CONXX+CONYY)+...
                            (CONXX+CONYY)*(PHIXX(L)+PHIYY(L)));
                        %}
                        SF(INOP(L))=SF(INOP(L))-W(J)*W(K)*DX*DY*...
                            (DCON*PHI(L)-DIFF*PHI(L)*(-1/(CON^2)+(1/10)*...
                            ((1-CON)^-2))*(CONX^2+CONY^2)-DIFF*PHI(L)*...
                            (1/CON+(1/10)*(1/(1-CON))-2*CHI)*...
                            (CONXX+CONYY)+(CONXX+CONYY)*...
                            (PHIXX(L)+PHIYY(L)));
                        for M = 1: 16
                            %{
                            SJ(INOP(L),INOP(M))=SJ(INOP(L),INOP(M))+W(J)...
                                *W(K)*DX*DY*(PHI(L)*PHI(M)/DT-DIFF...
                                *PHI(L)*(PHIXX(M)+PHIYY(M))+(PHIXX(L)+...
                                PHIYY(L))*(PHIXX(M)+PHIYY(M)));
                            %}
                            SJ(INOP(L),INOP(M))=SJ(INOP(L),INOP(M))+W(J)...
                                *W(K)*DX*DY*...
                                (PHI(L)*PHI(M)/DT-DIFF*PHI(L)*((2*PHI(M)...
                                /(CON^3)+2*PHI(M)/(10*(1-CON)^3))*...
                                (CONX^2+CONY^2)+(-1/CON^2+1/(10*(1-CON)^2))...
                                *(2*CONX*PHI(M)+2*CONY*PHI(M)))-DIFF*PHI(L)...
                                *((-PHI(M)/CON^2+PHI(M)/(10*(1-CON)^2))*...
                                (CONXX+CONYY)+(1/CON+1/(10*(1-CON))-2*CHI)...
                                *(PHIXX(M)+PHIYY(M)))+(PHIXX(M)+PHIYY(M))...
                                *(PHIXX(L)+PHIYY(L)));
                        end
                    end
                end
            end
        end
      
        display('SF SJ construction commences')
        %APPLY NATURAL BC==================================================
        display('Natural BC starts')
        [SF,SJ]=natural_BC(SF,SJ,NEX,NEY,NY,NFOUR);
        display('Natural BC commences')
        %INVERT JACOBIANS==================================================
        display('Sparsing')
        SJ=sparse(SJ);
        display('Matrix divisiion starts')
        sf=SJ\SF';
        SF=SJ\SF';
        SF=SF';
        display('Matrix divisiion commences')
        %THE SOLUTION VECTOR IS IN THE SF VECTOR NOW=======================
        JK = JK + 1;

        ERROR = sqrt(sum(SF.^2))
        
        C=C+SF;
        if ERROR<10^-6
            %C=C+SF;
            display('Newton-Raphson commences')
            if bypassed~=0
                bypassed=0;
            end
            break
        end
        
        if JK==1
            saved_error=ERROR;
        elseif JK>=20 || saved_error<=ERROR || isnan(ERROR)
            JK=100;
            display('BYPASSING!=========================================!')
            
            bypass_counter=bypass_counter+1
            
            bypassed=bypassed+1;
            
            break %bypass everything
            %error('Too many iterations')
        else
            saved_error=ERROR;
        end
    end
    
    %ADJUSTING TIME STEP TO MAKE CHANGE MORE APPARENT!=====================
    %[DTOO,DTO,DT]=update_time(DTO,DT,DT);
    if JK==100
        C=CO;
        CO=COO;
        next_step=DT/2;
        TIME=TIME-DTO;
        if bypassed==1
            DT=DTO;
            try
                DTO=DTOO;
            catch
                DTOO=DT;
            end
        elseif bypassed>=5
            display('NEWTON RAPHSON IS NOT CONVERGING! SIMULATION HAS BEEN DISRUPTED')
            break
            %error('Time step has been reduced too much!')
        end
        display('TIME STEP REDUCED BY ONE TENTH!!!')
        [DTOO,DTO,DT]=update_time(DTO,DT,next_step);
        TIME=TIME+DT
        continue
    elseif JK>=10
        display('HALVE TIME STEP!!!')
        [DTOO,DTO,DT]=update_time(DTO,DT,DT/2);
    elseif JK<=3
        display('INCREASE TIME STEP!!!')
        [DTOO,DTO,DT]=update_time(DTO,DT,DT*1.5);
        %{
        try
            abs_change=sum(abs(C-CO));
            if abs_change<standard_change
                [DTOO,DTO,DT]=update_time(DTO,DT,DT*1.5);
            else
                [DTOO,DTO,DT]=update_time(DTO,DT,DT);
            end
        catch
            standard_change=sum(abs(C-CO));
            [DTOO,DTO,DT]=update_time(DTO,DT,DT);
        end
        %}
    elseif DT<=DT_catchup
        display('CATCH UP TIME STEP!!!')
        %error('Just stop')
        [DTOO,DTO,DT]=update_time(DTO,DT,DT*catchup_pace);
    elseif max(abs(C-CO))<=0.02
        display('INCREASE TIME STEP!!!')
        %error('Just stop')
        [DTOO,DTO,DT]=update_time(DTO,DT,DT*1.5);
    else
        [DTOO,DTO,DT]=update_time(DTO,DT,DT);
    end
    %======================================================================
    
    %END OF NEWTON-RAPHSON=================================================
    
    %TAKE THE FRAME=======================================================
    solution=zeros(1,NX*NY);
    for i=1:NX*NY
        solution(i)=C((i-1)*4+1);
    end
    solution=reshape(solution,[NY,NX]);
    frame_index=frame_index+1;
    subplot(1,2,1)
    surf(x_coord,y_coord,solution);
    colorbar('southoutside')
    %caxis([min_bar max_bar]);
    axis([0 1 0 1 min_z max_z])
    xlabel('x'),ylabel('y'),zlabel('z')
    subplot(1,2,2)
    contourf(x_coord,y_coord,solution);
    %set(hC,'LineStyle','none');
    colorbar
    %caxis([min_bar max_bar]);
    suptitle('Solution to 2D nonlinear Cahn-Hilliard in unstable region')
    myClip(frame_index)=getframe(gcf);
    %======================================================================
    TIME=TIME+DTO
    if TIME>max_time
        break
    end
    %Measure time
    toc
end

bypass_counter
sf
display('Done!')

toc

%TAKE THE LAST FRAME=======================================================
%{
solution=zeros(1,NX*NY);
for i=1:NX*NY
    solution(i)=C((i-1)*4+1);
end
solution=reshape(solution,[NY,NX]);
frame_index=frame_index+1;
subplot(1,2,1)
surf(x_coord,y_coord,solution);
colorbar('southoutside')
%caxis([min_bar max_bar]);
axis([0 1 0 1 min_z max_z])
xlabel('x'),ylabel('y'),zlabel('z')
subplot(1,2,2)
contourf(x_coord,y_coord,solution);
%set(hC,'LineStyle','none');
colorbar
%caxis([min_bar max_bar]);
suptitle('Solution to 2D nonlinear Cahn-Hilliard in unstable region')
myClip(frame_index)=getframe(gcf);
%}
%==========================================================================

figure(1005)
surf(x_coord,y_coord,solution)

final_fig=figure(200);
set(final_fig, 'Position', pos_size)

movie(gcf,myClip,5)

%%{
%EXPORT THE MOVIE================================
%movie2avi(myClip,'Solution to 2D nonlinear Cahn-Hilliard in unstable region (Fixed axis)')
video=VideoWriter('10201 nodes.avi', 'Uncompressed AVI');
video.FrameRate=12;
open(video)
writeVideo(video,myClip);
close(video)
%%}
%================================

toc

end

%==========================================================================
%==========================================================================
%==========================================================================

function [DTOO,DTO,DT]=update_time(DTO,DT,next_time_step)

DTOO=DTO;
DTO=DT;
DT=next_time_step;
end

%==========================================================================
%==========================================================================
%==========================================================================

function [SF,SJ]=natural_BC(SF,SJ,NEX,NEY,NY,NFOUR)
% APPLY NATURAL BOUNDARY CONDITIONS

display('Loop 1')
%Left
for J = 2:4: 2 + 4 * NEY
    %J1 = J + 1;
    %J2 = J + 2;
    SJ(J, :) = zeros(1,NFOUR);
    %SJ(J1, :) = zeros(1,NFOUR);
    %SJ(J2, :) = zeros(1,NFOUR);
    SF(J) = 0;
    %SF(J1) = 0;
    %SF(J2) = 0;
    SJ(J, J) = 1;
    %SJ(J1, J1) = 1;
    %SJ(J2, J2) = 1;
end
display('Loop 2')
%Bottom
for J = 2:4 * NY: 2 + 4 * NEX * NY
    %display('stage 1')
    J1 = J + 1;
    %J2 = J + 2;
    %display('stage 2')
    %SJ(J, :) = zeros(1,NFOUR);
    SJ(J1, :) = zeros(1,NFOUR);
    %SJ(J2, :) = zeros(1,NFOUR);
    %display('stage 3')
    %SF(J) = 0;
    SF(J1) = 0;
    %SF(J2) = 0;
    %SJ(J, J) = 1;
    SJ(J1, J1) = 1;
    %SJ(J2, J2) = 1;
end
display('Loop 3')
%Top
for J = 2 + 4 * NEY:4 * NY: 2 + 4 * NEY + 4 * NEX * NY
    J1 = J + 1;
    %J2 = J + 2;
    %SJ(J, :) = zeros(1,NFOUR);
    SJ(J1, :) = zeros(1,NFOUR);
    %SJ(J2, :) = zeros(1,NFOUR);
    %SF(J) = 0;
    SF(J1) = 0;
    %SF(J2) = 0;
    %=======================================================================6
    %SJ(J, J) = 1.000;
    SJ(J1,J1)=1;
    %SJ(J2, J2)=1;
end
display('Loop 4')
%Right
for J = 2 + 4*NEX*NY:4:2+4*NEY+4*NEX*NY
    %J1=J+1;
    %J2 = J + 2;
    SJ(J, :) = zeros(1,NFOUR);
    %SJ(J1, :) = zeros(1,NFOUR);
    %SJ(J2, :) = zeros(1,NFOUR);
    SF(J) = 0;
    %SF(J1) = 0;
    %SF(J2) = 0;
    SJ(J, J) = 1;
    %SJ(J1, J1) = 1;
    %SJ(J2, J2) = 1;
end
end

%==========================================================================
%==========================================================================
%==========================================================================

function solution=F1(Z)
solution=1-3*Z^2+2*Z^3;
end
function solution=F2(Z,DZ)
solution=(Z-2*Z^2+Z^3)*DZ;
end
function solution=F3(Z)
solution=3*Z^2-2*Z^3;
end
function solution=F4(Z,DZ)
solution=(Z^3-Z^2)*DZ;
end
function solution=F1Z(Z,DZ)
solution=(-6*Z+6*Z^2)/DZ;
end
function solution=F2Z(Z)
solution=1-4*Z+3*Z^2;
end
function solution=F3Z(Z,DZ)
solution=(6*Z-6*Z^2)/DZ;
end
function solution=F4Z(Z)
solution=3*Z^2-2*Z;
end
function solution=F1ZZ(Z,DZ)
solution=(-6+12*Z)/DZ^2;
end
function solution=F2ZZ(Z,DZ)
solution=(-4+6*Z)/DZ;
end
function solution=F3ZZ(Z,DZ)
solution=(6-12*Z)/DZ^2;
end
function solution=F4ZZ(Z,DZ)
solution=(6*Z-2)/DZ;
end

function solution=fPHI(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1(GP1)*F1(GP2);
solution(2)=F2(GP1,DX)*F1(GP2);
solution(3)=F1(GP1)*F2(GP2,DY);
solution(4)=F2(GP1,DX)*F2(GP2,DY);
solution(5)=F1(GP1)*F3(GP2);
solution(6)=F2(GP1,DX)*F3(GP2);
solution(7)=F1(GP1)*F4(GP2,DY);
solution(8)=F2(GP1,DX)*F4(GP2,DY);
solution(9)=F3(GP1)*F1(GP2);
solution(10)=F4(GP1,DX)*F1(GP2);
solution(11)=F3(GP1)*F2(GP2,DY);
solution(12)=F4(GP1,DX)*F2(GP2,DY);
solution(13)=F3(GP1)*F3(GP2);
solution(14)=F4(GP1,DX)*F3(GP2);
solution(15)=F3(GP1)*F4(GP2,DY);
solution(16)=F4(GP1,DX)*F4(GP2,DY);
end

function solution=fPHIX(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1Z(GP1,DX)*F1(GP2);
solution(2)=F2Z(GP1)*F1(GP2);
solution(3)=F1Z(GP1,DX)*F2(GP2,DY);
solution(4)=F2Z(GP1)*F2(GP2,DY);
solution(5)=F1Z(GP1,DX)*F3(GP2);
solution(6)=F2Z(GP1)*F3(GP2);
solution(7)=F1Z(GP1,DX)*F4(GP2,DY);
solution(8)=F2Z(GP1)*F4(GP2,DY);
solution(9)=F3Z(GP1,DX)*F1(GP2);
solution(10)=F4Z(GP1)*F1(GP2);
solution(11)=F3Z(GP1,DX)*F2(GP2,DY);
solution(12)=F4Z(GP1)*F2(GP2,DY);
solution(13)=F3Z(GP1,DX)*F3(GP2);
solution(14)=F4Z(GP1)*F3(GP2);
solution(15)=F3Z(GP1,DX)*F4(GP2,DY);
solution(16)=F4Z(GP1)*F4(GP2,DY);
end

function solution=fPHIY(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1(GP1)*F1Z(GP2,DY);
solution(2)=F2(GP1,DX)*F1Z(GP2,DY);
solution(3)=F1(GP1)*F2Z(GP2);
solution(4)=F2(GP1,DX)*F2Z(GP2);
solution(5) =F1(GP1)*F3Z(GP2,DY);
solution(6)=F2(GP1,DX)*F3Z(GP2,DY);
solution(7)=F1(GP1)*F4Z(GP2);
solution(8)=F2(GP1,DX)*F4Z(GP2);
solution(9)=F3(GP1)*F1Z(GP2,DY);
solution(10)=F4(GP1,DX)*F1Z(GP2,DY);
solution(11)=F3(GP1)*F2Z(GP2);
solution(12)=F4(GP1,DX)*F2Z(GP2);
solution(13)=F3(GP1)*F3Z(GP2,DY);
solution(14)=F4(GP1,DX)*F3Z(GP2,DY);
solution(15)=F3(GP1)*F4Z(GP2);
solution(16)=F4(GP1,DX)*F4Z(GP2);
end

function solution=fPHIXY(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1Z(GP1,DX)*F1Z(GP2,DY);
solution(2)=F2Z(GP1)*F1Z(GP2,DY);
solution(3)=F1Z(GP1,DX)*F2Z(GP2);
solution(4)=F2Z(GP1)*F2Z(GP2);
solution(5)=F1Z(GP1,DX)*F3Z(GP2,DY);
solution(6)=F2Z(GP1)*F3Z(GP2,DY);
solution(7)=F1Z(GP1,DX)*F4Z(GP2);
solution(8)=F2Z(GP1)*F4Z(GP2);
solution(9)=F3Z(GP1,DX)*F1Z(GP2,DY);
solution(10)=F4Z(GP1)*F1Z(GP2,DY);
solution(11) = F3Z(GP1,DX)*F2Z(GP2);
solution(12)=F4Z(GP1)*F2Z(GP2);
solution(13)=F3Z(GP1,DX)*F3Z(GP2,DY);
solution(14)=F4Z(GP1)*F3Z(GP2,DY);
solution(15)=F3Z(GP1,DX)*F4Z(GP2);
solution(16)=F4Z(GP1)*F4Z(GP2);
end

function solution=fPHIXX(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1ZZ(GP1,DX)*F1(GP2);
solution(2)=F2ZZ(GP1,DX)*F1(GP2);
solution(3)=F1ZZ(GP1,DX)*F2(GP2,DY);
solution(4)=F2ZZ(GP1,DX)*F2(GP2,DY);
solution(5)=F1ZZ(GP1,DX)*F3(GP2);
solution(6)=F2ZZ(GP1,DX)*F3(GP2);
solution(7)=F1ZZ(GP1,DX)*F4(GP2,DY);
solution(8)=F2ZZ(GP1,DX)*F4(GP2,DY);
solution(9)=F3ZZ(GP1,DX)*F1(GP2);
solution(10)=F4ZZ(GP1,DX)*F1(GP2);
solution(11)=F3ZZ(GP1,DX)*F2(GP2,DY);
solution(12)=F4ZZ(GP1,DX)*F2(GP2,DY);
solution(13)=F3ZZ(GP1,DX)*F3(GP2);
solution(14)=F4ZZ(GP1,DX)*F3(GP2);
solution(15)=F3ZZ(GP1,DX)*F4(GP2,DY);
solution(16)=F4ZZ(GP1,DX)*F4(GP2,DY);
end

function solution=fPHIYY(GP1,GP2,DX,DY)
solution=zeros(1,16);
solution(1)=F1(GP1)*F1ZZ(GP2,DY);
solution(2)=F2(GP1,DX)*F1ZZ(GP2,DY);
solution(3)=F1(GP1)*F2ZZ(GP2,DY);
solution(4)=F2(GP1,DX)*F2ZZ(GP2,DY);
solution(5)=F1(GP1)*F3ZZ(GP2,DY);
solution(6)=F2(GP1,DX)*F3ZZ(GP2,DY);
solution(7)=F1(GP1)*F4ZZ(GP2,DY);
solution(8)=F2(GP1,DX)*F4ZZ(GP2,DY);
solution(9)=F3(GP1)*F1ZZ(GP2,DY);
solution(10)=F4(GP1,DX)*F1ZZ(GP2,DY);
solution(11)=F3(GP1)*F2ZZ(GP2,DY);
solution(12)=F4(GP1,DX)*F2ZZ(GP2,DY);
solution(13)=F3(GP1)*F3ZZ(GP2,DY);
solution(14)=F4(GP1,DX)*F3ZZ(GP2,DY);
solution(15)=F3(GP1)*F4ZZ(GP2,DY);
solution(16)=F4(GP1,DX)*F4ZZ(GP2,DY);
end

function [PHI,PHIX,PHIY,PHIXX,PHIYY,PHIXY]=TFUNCT(GP1,GP2,DX,DY)
PHI=fPHI(GP1,GP2,DX,DY);
PHIX=fPHIX(GP1,GP2,DX,DY);
PHIY=fPHIY(GP1,GP2,DX,DY);
PHIXX=fPHIXX(GP1,GP2,DX,DY);
PHIYY=fPHIYY(GP1,GP2,DX,DY);
PHIXY=fPHIXY(GP1,GP2,DX,DY);
end