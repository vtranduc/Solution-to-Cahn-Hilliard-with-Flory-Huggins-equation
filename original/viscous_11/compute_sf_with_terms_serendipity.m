function sf=compute_sf_with_terms_serendipity(nfour,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,determinants,diffT,wTerms,terms)

%TO BE FILLED IN LATER!

end

function solution=compute_sfth_sf_with_terms_serendipity(sfth,conc_,conco_,nx,ny,nz,nex,ney,weights,dt,diffT,determinants,wTerms,terms)

solution=0;

w=[5/18 4/9 5/18];

[node,type]=analyze_gbs_3d(sfth);
[xth,yth,zth]=get_xyzth_3d(node,nx,ny);
elements=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);

for i_positional_element=1:1:8
    if elements(i_positional_element)~=0
        
        con=dim_reducer2_3d(conc_(elements(i_positional_element),1,:,:,:));
        cono=dim_reducer1_3d(conco_(elements(i_positional_element),:,:,:));
        orientation=9-i_positional_element;
        cont=(con-cono)/dt;
        
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    solution=solution+w(ix)*w(iy)*w(iz)...
                        *determinants(elements(i_positional_element),ix,iy,iz)*(...
                        weights(elements(i_positional_element),orientation,type,1,ix,iy,iz)*cont(ix,iy,iz)...
                        -diffT*weights(elements(i_positional_element),orientation,type,1,ix,iy,iz)*...
                        (terms(elements(i_positional_element),1,ix,iy,iz)+...
                        terms(elements(i_positional_element),2,ix,iy,iz))+...
                        terms(elements(i_positional_element),3,ix,iy,iz)*...
                        wTerms(elements(i_positional_element),orientation,type,ix,iy,iz)...
                        );
                end
            end
        end
    end
end
% solution=dxyz*solution;
end