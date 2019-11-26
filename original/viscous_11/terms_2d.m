function [terms_sf,terms_sj]=terms_2d(c,ne,n1,n2,Two_chi_n1,...
    grad_T,coef_T,diff,terms_assist,thermophoresis,...
    entropy)

if grad_T==0
    %=====================MUST FIX THIS==============================
%     terms=zeros(ne,3,3,6);

    terms_sf=zeros(ne,3,3,2);
    terms_sj=zeros(ne,3,3,6);
    
    for e=1:1:ne
%         [inx,~]=inxiny_elemental(e,ney);
        for ix=1:1:3
            for iy=1:1:3
                
                sub=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                c2=1-c(ix,iy,1,e);
                
                terms_sf(e,ix,iy,2)=c(ix,iy,4,e)+c(ix,iy,5,e);
                
                terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1;
                terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                terms_sj(e,ix,iy,3)=terms_sf(e,ix,iy,2);
                terms_sj(e,ix,iy,4)=...
                    2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)*sub;
                terms_sj(e,ix,iy,5)=c(ix,iy,2,e);
                terms_sj(e,ix,iy,6)=c(ix,iy,3,e);
                
                terms_sf(e,ix,iy,1)=coef_T*(terms_sj(e,ix,iy,1)*terms_sj(e,ix,iy,3)...
                    +terms_sj(e,ix,iy,2)*sub);
                
%                 sub1=-(c(ix,iy,1,e)^2*n1)^-1+((1-c(ix,iy,1,e))^2*n2)^-1;
%                 sub2=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
%                 sub3=(c(ix,iy,1,e)*n1)^-1+((1-c(ix,iy,1,e))*n2)^-1-Two_chi_n1;
%                 sub4=c(ix,iy,4,e)+c(ix,iy,5,e);
%                 terms(e,ix,iy,1)=sub1*sub2;
%                 terms(e,ix,iy,2)=sub3*sub4;
%                 terms(e,ix,iy,3)=sub4;
%                 terms(e,ix,iy,4)=...
%                     2*((c(ix,iy,1,e)^3*n1)^-1+((1-c(ix,iy,1,e))^3*n2)^-1)...
%                     *sub2;
%                 terms(e,ix,iy,5)=sub1;
%                 terms(e,ix,iy,6)=sub3;
            end
        end
    end
    %==================================================================
elseif grad_T==1

    if thermophoresis==0
        
%         error('fdadfa')
        
        terms_sf=zeros(ne,3,3,4);
        terms_sj=zeros(ne,3,3,8);
    
        for e=1:1:ne
            
            for ix=1:1:3
                for iy=1:1:3

                    sub1=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                    c2=1-c(ix,iy,1,e);

                    sub2=diff*(...
                        (log(c(ix,iy,1,e))+1)/n1...
                        -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                        *(0.5-c(ix,iy,1,e))...
                        +terms_assist(e,ix,iy)*(c(ix,iy,1,e)-0.5)...
                        );
                    
                    
                    %-----------------------------
%                     a=(...
%                         (log(c(ix,iy,1,e))+1)/n1...
%                         -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
%                         *(0.5-c(ix,iy,1,e)));
%                     b=terms_assist(e,ix,iy)*(c(ix,iy,1,e)-0.5);
%                     [a b]
                    %------------------------------
                    
                    
                    terms_sf(e,ix,iy,2)=sub2*coef_T(e,ix,iy,2);
                    terms_sf(e,ix,iy,3)=sub2*coef_T(e,ix,iy,3);
                    terms_sf(e,ix,iy,4)=c(ix,iy,4,e)+c(ix,iy,5,e);
                    

                    terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy);
                    terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                    terms_sj(e,ix,iy,3)=(coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                        +coef_T(e,ix,iy,3)*c(ix,iy,2,3))...
                        +coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4);

                    terms_sj(e,ix,iy,4)=diff...
                        *2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)...
                        *coef_T(e,ix,iy,1)*sub1;

                    sub3=diff*(terms_sj(e,ix,iy,1)+terms_assist(e,ix,iy));

                    terms_sj(e,ix,iy,5)=sub3*coef_T(e,ix,iy,2);
                    terms_sj(e,ix,iy,6)=sub3*coef_T(e,ix,iy,3);
                    sub4=2*coef_T(e,ix,iy,1);
                    terms_sj(e,ix,iy,7)=sub4*c(ix,iy,2,e);
                    terms_sj(e,ix,iy,8)=sub4*c(ix,iy,3,e);

                    dot_Tc=coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                        +coef_T(e,ix,iy,3)*c(ix,iy,2,3);

                    terms_sf(e,ix,iy,1)=diff*(terms_sj(e,ix,iy,1)*(...
                        dot_Tc+coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4))...
                        ...
                        +terms_sj(e,ix,iy,2)*coef_T(e,ix,iy,1)*sub1...
                        +terms_assist(e,ix,iy)*dot_Tc);

                    terms_sj(e,ix,iy,1)=diff*terms_sj(e,ix,iy,1)*coef_T(e,ix,iy,1);
                    terms_sj(e,ix,iy,2)=diff*terms_sj(e,ix,iy,2); 


                end
            end
        end
        
    elseif thermophoresis~=0
        
%         error('toch')
        
        terms_sf=zeros(ne,3,3,4);
        terms_sj=zeros(ne,3,3,10);
        
        for e=1:1:ne
    %         [inx,~]=inxiny_elemental(e,ney);
            for ix=1:1:3
                for iy=1:1:3

                    sub1=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                    c2=1-c(ix,iy,1,e);

                    sub2=diff*(...
                        (log(c(ix,iy,1,e))+1)/n1...
                        -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                        *(0.5-c(ix,iy,1,e))...
                        +terms_assist(e,ix,iy)*(c(ix,iy,1,e)-0.5)...
                        ...
                        ...%==============THERMOPHORESIS============
                        )+thermophoresis*c(ix,iy,1,e)*c2...
                        ...%=========================================
                        ...
                        ;
                    
                    
                    terms_sf(e,ix,iy,2)=sub2*coef_T(e,ix,iy,2);
                    terms_sf(e,ix,iy,3)=sub2*coef_T(e,ix,iy,3);
                    terms_sf(e,ix,iy,4)=c(ix,iy,4,e)+c(ix,iy,5,e);

                    terms_sj(e,ix,iy,1)=(c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy);
                    terms_sj(e,ix,iy,2)=-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1;
                    terms_sj(e,ix,iy,3)=(coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                        +coef_T(e,ix,iy,3)*c(ix,iy,2,3))...
                        +coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4);

                    terms_sj(e,ix,iy,4)=diff...
                        *2*((c(ix,iy,1,e)^3*n1)^-1+(c2^3*n2)^-1)...
                        *coef_T(e,ix,iy,1)*sub1;

                    sub3=diff*(terms_sj(e,ix,iy,1)+terms_assist(e,ix,iy));

                    terms_sj(e,ix,iy,5)=sub3*coef_T(e,ix,iy,2);
                    terms_sj(e,ix,iy,6)=sub3*coef_T(e,ix,iy,3);
                    sub4=2*coef_T(e,ix,iy,1);
                    terms_sj(e,ix,iy,7)=sub4*c(ix,iy,2,e);
                    terms_sj(e,ix,iy,8)=sub4*c(ix,iy,3,e);

                    dot_Tc=coef_T(e,ix,iy,2)*c(ix,iy,2,e)...
                        +coef_T(e,ix,iy,3)*c(ix,iy,2,3);

                    terms_sf(e,ix,iy,1)=diff*(terms_sj(e,ix,iy,1)*(...
                        dot_Tc+coef_T(e,ix,iy,1)*terms_sf(e,ix,iy,4))...
                        ...
                        +terms_sj(e,ix,iy,2)*coef_T(e,ix,iy,1)*sub1...
                        +terms_assist(e,ix,iy)*dot_Tc);

                    terms_sj(e,ix,iy,1)=diff*terms_sj(e,ix,iy,1)*coef_T(e,ix,iy,1);
                    terms_sj(e,ix,iy,2)=diff*terms_sj(e,ix,iy,2);
                    
                    
                    %=============THERMOPHORESIS================
                    sub5=sub3+thermophoresis*(1-2*c(ix,iy,1,e));
                    
                    terms_sj(e,ix,iy,9)=sub5*coef_T(e,ix,iy,2);
                    terms_sj(e,ix,iy,10)=sub5*coef_T(e,ix,iy,3);
                    %===========================================


                end
            end
        end
        
    
    end
    
elseif grad_T==2
    
%     error('dfafd')
    
    terms_sf=zeros(ne,3,3,5);
    terms_sj=zeros(ne,3,3,9);
    
    terms_assist_=diff*terms_assist;

    for e=1:1:ne
%         [inx,~]=inxiny_elemental(e,ney);
        for ix=1:1:3
            for iy=1:1:3

%                 sub1=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;

                c2=1-c(ix,iy,1,e);
                c1c2=c(ix,iy,1,e)*c2;
                
                sub1=1-2*c(ix,iy,1,e);
                
                sub2=diff...
                    *((log(c(ix,iy,1,e))+1)/n1...
                    -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                    *(0.5-c(ix,iy,1,e))...
                    -entropy*sub1/(n1*coef_T(e,ix,iy,1)));
                
                sub3=thermophoresis*c1c2;
                
                sub4=c1c2*sub2+sub3;
                
                terms_sf(e,ix,iy,1)=sub4*coef_T(e,ix,iy,2);
                terms_sf(e,ix,iy,2)=sub4*coef_T(e,ix,iy,3);
                
                sub5=diff...
                    *((c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy));
                
                laplacian=c(ix,iy,4,e)+c(ix,iy,5,e);

                sub6=c1c2*coef_T(e,ix,iy,1)*sub5+laplacian*sub1;
                
                terms_sf(e,ix,iy,3)=sub6*c(ix,iy,2,e);
                terms_sf(e,ix,iy,4)=sub6*c(ix,iy,3,e);
                
                terms_sf(e,ix,iy,5)=c1c2*laplacian;
                
                
                %----------------------------------------------------
                
                sub7=sub1*sub2...
                    +c1c2*(sub5+terms_assist_(e,ix,iy))...
                    +thermophoresis*sub1;
                
                terms_sj(e,ix,iy,1)=sub7*coef_T(e,ix,iy,2);
                terms_sj(e,ix,iy,2)=sub7*coef_T(e,ix,iy,3);
                
                terms_sj(e,ix,iy,3)=...
                    coef_T(e,ix,iy,1)*c1c2*sub5...
                    +laplacian*sub1;
                
                
                sub8=...
                    coef_T(e,ix,iy,1)...
                    *(sub1*sub5...
                    +diff*c1c2*(-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1))...
                    -2*laplacian;
                
                terms_sj(e,ix,iy,4)=sub8*c(ix,iy,2,e);
                terms_sj(e,ix,iy,5)=sub8*c(ix,iy,3,e);
                
                terms_sj(e,ix,iy,6)=sub1*c(ix,iy,2,e);
                terms_sj(e,ix,iy,7)=sub1*c(ix,iy,3,e);
                
                terms_sj(e,ix,iy,8)=sub1*laplacian;
                terms_sj(e,ix,iy,9)=c1c2;


            end
        end
    end
    
elseif grad_T==3
    
    terms_sf=zeros(ne,3,3,5);
    terms_sj=zeros(ne,3,3,13);
    
    terms_assist_=diff*terms_assist;
    
    %---------
    significance=0.0001;
    %-----------

    for e=1:1:ne
%         [inx,~]=inxiny_elemental(e,ney);
        for ix=1:1:3
            for iy=1:1:3
                
                c2=1-c(ix,iy,1,e);
                c1c2=c(ix,iy,1,e)*c2;
                laplacian=c(ix,iy,4,e)+c(ix,iy,5,e);
                
                terms_sj(e,ix,iy,12)=1/(c(ix,iy,1,e)*c2);
                
                %--------------
%                 terms_sj(e,ix,iy,12)=1;
                terms_sj(e,ix,iy,12)=significance*terms_sj(e,ix,iy,12)+1;
                %----------------
                
                sub1=1-2*c(ix,iy,1,e);
                sub2=diff...
                    *((log(c(ix,iy,1,e))+1)/n1...
                    -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
                    *(0.5-c(ix,iy,1,e))...
                    -entropy*sub1/(n1*coef_T(e,ix,iy,1)))...
                    +thermophoresis*c1c2;
                
                terms_sf(e,ix,iy,1)=sub2*coef_T(e,ix,iy,2);
                terms_sf(e,ix,iy,2)=sub2*coef_T(e,ix,iy,3);
                
                dk_dc=(2*c(ix,iy,1,e)-1)/(c(ix,iy,1,e)*c2)^2;
                
                d2k_dc2=2*(((c(ix,iy,1,e)^3)*c2)^-1-...
                    (c(ix,iy,1,e)*c2)^-2+...
                    (c(ix,iy,1,e)*c2^3)^-1);
                
                %-----------------
%                 dk_dc=0;
%                 d2k_dc2=0;

                dk_dc=significance*dk_dc;
                d2k_dc2=significance*d2k_dc2;
                %--------------------
                
                sub3=diff...
                    *((c(ix,iy,1,e)*n1)^-1+...
                    (c2*n2)^-1-Two_chi_n1(e,ix,iy));
                
                grad_squared=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
                
                sub4=coef_T(e,ix,iy,1)*sub3...
                    -0.5*d2k_dc2*grad_squared...
                    -dk_dc*laplacian;
                
                terms_sf(e,ix,iy,3)=sub4*c(ix,iy,2,e);
                terms_sf(e,ix,iy,4)=sub4*c(ix,iy,3,e);
                
                terms_sf(e,ix,iy,5)=terms_sj(e,ix,iy,12)*laplacian;
                
                
                
                sub5=sub3+terms_assist_(e,ix,iy)...
                    +thermophoresis*sub1;
                
                terms_sj(e,ix,iy,1)=sub5*coef_T(e,ix,iy,2);
                terms_sj(e,ix,iy,2)=sub5*coef_T(e,ix,iy,3);
                
                terms_sj(e,ix,iy,3)=sub4;
                
                d3k_dc3_two=3*(-(c(ix,iy,1,e)^4*c2)^-1+...
                    (c(ix,iy,1,e)^3*c2^2)^-1-...
                    (c(ix,iy,1,e)^2*c2^3)^-1+...
                    (c(ix,iy,1,e)*c2^4)^-1);
                
                %-------------
%                 d3k_dc3_two=0;

                d3k_dc3_two=significance*d3k_dc3_two;
                %------------
                
                sub6=diff*coef_T(e,ix,iy,1)...
                    *(-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1)...
                    -d3k_dc3_two*grad_squared...
                    -d2k_dc2*laplacian;
                
                
                terms_sj(e,ix,iy,4)=sub6*c(ix,iy,2,e);
                terms_sj(e,ix,iy,5)=sub6*c(ix,iy,3,e);
                
                terms_sj(e,ix,iy,6)=d2k_dc2*c(ix,iy,2,e);
                terms_sj(e,ix,iy,7)=d2k_dc2*c(ix,iy,3,e);
                
                terms_sj(e,ix,iy,8)=dk_dc*c(ix,iy,2,e);
                terms_sj(e,ix,iy,9)=dk_dc*c(ix,iy,3,e);
                
                terms_sj(e,ix,iy,10)=c(ix,iy,2,e);
                terms_sj(e,ix,iy,11)=c(ix,iy,3,e);
                
                terms_sj(e,ix,iy,13)=dk_dc*laplacian;
                

% %                 sub1=c(ix,iy,2,e)^2+c(ix,iy,3,e)^2;
% 
%                 c2=1-c(ix,iy,1,e);
%                 c1c2=c(ix,iy,1,e)*c2;
%                 
%                 sub1=1-2*c(ix,iy,1,e);
%                 
%                 sub2=diff...
%                     *((log(c(ix,iy,1,e))+1)/n1...
%                     -(log(c2)+1)/n2+Two_chi_n1(e,ix,iy)...
%                     *(0.5-c(ix,iy,1,e))...
%                     -entropy*sub1/(n1*coef_T(e,ix,iy,1)));
%                 
%                 sub3=thermophoresis*c1c2;
%                 
%                 sub4=c1c2*sub2+sub3;
%                 
%                 terms_sf(e,ix,iy,1)=sub4*coef_T(e,ix,iy,2);
%                 terms_sf(e,ix,iy,2)=sub4*coef_T(e,ix,iy,3);
%                 
%                 sub5=diff...
%                     *((c(ix,iy,1,e)*n1)^-1+(c2*n2)^-1-Two_chi_n1(e,ix,iy));
%                 
%                 laplacian=c(ix,iy,4,e)+c(ix,iy,5,e);
% 
%                 sub6=c1c2*coef_T(e,ix,iy,1)*sub5+laplacian*sub1;
%                 
%                 terms_sf(e,ix,iy,3)=sub6*c(ix,iy,2,e);
%                 terms_sf(e,ix,iy,4)=sub6*c(ix,iy,3,e);
%                 
%                 terms_sf(e,ix,iy,5)=c1c2*laplacian;
%                 
%                 
%                 %----------------------------------------------------
%                 
%                 sub7=sub1*sub2...
%                     +c1c2*(sub5+terms_assist_(e,ix,iy))...
%                     +thermophoresis*sub1;
%                 
%                 terms_sj(e,ix,iy,1)=sub7*coef_T(e,ix,iy,2);
%                 terms_sj(e,ix,iy,2)=sub7*coef_T(e,ix,iy,3);
%                 
%                 terms_sj(e,ix,iy,3)=...
%                     coef_T(e,ix,iy,1)*c1c2*sub5...
%                     +laplacian*sub1;
%                 
%                 
%                 sub8=...
%                     coef_T(e,ix,iy,1)...
%                     *(sub1*sub5...
%                     +diff*c1c2*(-(c(ix,iy,1,e)^2*n1)^-1+(c2^2*n2)^-1))...
%                     -2*laplacian;
%                 
%                 terms_sj(e,ix,iy,4)=sub8*c(ix,iy,2,e);
%                 terms_sj(e,ix,iy,5)=sub8*c(ix,iy,3,e);
%                 
%                 terms_sj(e,ix,iy,6)=sub1*c(ix,iy,2,e);
%                 terms_sj(e,ix,iy,7)=sub1*c(ix,iy,3,e);
%                 
%                 terms_sj(e,ix,iy,8)=sub1*laplacian;
%                 terms_sj(e,ix,iy,9)=c1c2;


            end
        end
    end


end
    
end