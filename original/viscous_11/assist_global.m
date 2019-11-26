function terms_assist=assist_global(...
    grad_T,ne,temperature,entropy,n1)

if grad_T==0
    terms_assist=NaN;
%     sj_assist=NaN;
elseif grad_T==1
    terms_assist=zeros(ne,3,3);
%     sj_assist=zeros(ne,16,3,3,2);

    c1=2*entropy/n1;

    for e=1:1:ne
        for ix=1:1:3
            for iy=1:1:3
                terms_assist(e,ix,iy)=c1/temperature(e,ix,iy,1);

%                 for jlocal=1:1:16
% 
%                     sj_assist(e,jlocal,ix,iy,1)=...
%                       temperature(e,ix,iy,2)*weights(ix,iy,jlocal,2)...
%                       +temperature(e,ix,iy,3)*weights(ix,iy,jlocal,3)...
%                       +temperature(e,ix,iy,1)*wTerms(ix,iy,jlocal);
% 
% 
%                     sj_assist(e,jlocal,ix,iy,2)=...
%                         terms_assist(e,ix,iy)...
%                         *temperature(e,ix,iy,2)*weights(ix,iy,jlocal,2)...
%                         +terms_assist(e,ix,iy)...
%                         *temperature(e,ix,iy,3)*weights(ix,iy,jlocal,3);
%                 end

            end
        end
    end
end

end