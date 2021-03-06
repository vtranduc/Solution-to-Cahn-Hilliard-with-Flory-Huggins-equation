function sj=test26(nnz_,dt,weights,dxyz,diffT,wTerms,terms,sj_properties)

sj=zeros(1,nnz_);

w=[5/18 4/9 5/18];

if isscalar(diffT)
    parfor sjth=1:1:nnz_
%         sj(sjth)=compute_sjth_sj_with_terms_and_sj_properties_3d(sjth,...
%             w,dt,diffT,weights,dxyz,wTerms,terms,sj_properties);
    end
else
    parfor sjth=1:1:nnz_
        sjth_properties=sj_properties(sjth,:,:);
        sj(sjth)=...
            warui(...
            w,dt,diffT,weights,dxyz,wTerms,terms,sjth_properties);
    end
end

end

% === Temperature gradient is present =============================================

function solution=warui(w,dt,diffT,weights,dxyz,wTerms,terms,sjth_properties)
% disp('hefda')
solution=0;
for i_positional_element=1:1:8
    if sjth_properties(1,i_positional_element,1)~=0
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        weights(sjth_properties(1,i_positional_element,2),sjth_properties(1,i_positional_element,4),1,ix,iy,iz)*weights(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),1,ix,iy,iz)/dt...
                        -diffT(sjth_properties(1,i_positional_element,6),ix)*weights(sjth_properties(1,i_positional_element,2),sjth_properties(1,i_positional_element,4),1,ix,iy,iz)...
                        *(weights(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),1,ix,iy,iz)...
                        *terms(sjth_properties(1,i_positional_element,1),4,ix,iy,iz)...
                        +weights(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),1,ix,iy,iz)...
                        *terms(sjth_properties(1,i_positional_element,1),5,ix,iy,iz)...
                        +terms(sjth_properties(1,i_positional_element,1),6,ix,iy,iz)...
                        *wTerms(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),ix,iy,iz))...
                        +wTerms(sjth_properties(1,i_positional_element,2),sjth_properties(1,i_positional_element,4),ix,iy,iz)...
                        *wTerms(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),ix,iy,iz)...
                        ...
                        );
                end
            end
        end
    end
end
solution=solution*dxyz;
end

% === No temperature gradient =============================================

% function solution=compute_sjth_sj_with_terms_and_sj_properties_3d(...
%     sjth,w,dt,diffT,weights,dxyz,wTerms,terms,sj_properties)
% 
% solution=0;
% for i_positional_element=1:1:8
%     if sjth_properties(1,i_positional_element,1)~=0
%         for ix=1:1:3
%             for iy=1:1:3
%                 for iz=1:1:3
%                     
%                     solution=solution+w(ix)*w(iy)*w(iz)*(...
%                         weights(sjth_properties(1,i_positional_element,2),sjth_properties(1,i_positional_element,4),1,ix,iy,iz)*weights(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),1,ix,iy,iz)/dt...
%                         -diffT*weights(sjth_properties(1,i_positional_element,2),sjth_properties(1,i_positional_element,4),1,ix,iy,iz)...
%                         *(weights(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),1,ix,iy,iz)...
%                         *terms(sjth_properties(1,i_positional_element,1),4,ix,iy,iz)...
%                         +weights(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),1,ix,iy,iz)...
%                         *terms(sjth_properties(1,i_positional_element,1),5,ix,iy,iz)...
%                         +terms(sjth_properties(1,i_positional_element,1),6,ix,iy,iz)...
%                         *wTerms(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),ix,iy,iz))...
%                         +wTerms(sjth_properties(1,i_positional_element,2),sjth_properties(1,i_positional_element,4),ix,iy,iz)...
%                         *wTerms(sjth_properties(1,i_positional_element,3),sjth_properties(1,i_positional_element,5),ix,iy,iz)...
%                         ...
%                         );
%                     
%                 end
%             end
%         end
%     end
% end
% solution=solution*dxyz;
% end