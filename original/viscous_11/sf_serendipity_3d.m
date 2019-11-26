function sf=sf_serendipity_3d(...
    nfour,nx,ny,nz,nex,ney,weights,dxyz,...
    wTerms,terms)

sf=zeros(1,nfour);

w=[5/18 4/9 5/18];

parfor gbf=1:1:nfour %PARFOR MUST BE FEASIBLE!!!
    sf(gbf)=compute_sfth_sf_with_terms_serendipity(...
        gbf,nx,ny,nz,nex,ney,weights,dxyz,wTerms,terms,w);
end

end

function solution=compute_sfth_sf_with_terms_serendipity(...
    gbf,nx,ny,nz,nex,ney,weights,dxyz,wTerms,terms,w)

solution=0;

[node,type]=analyze_gbs_serendipity(gbf);

[xth,yth,zth]=get_xyzth_3d(node,nx,ny);
elements=get_elements_3d(xth,yth,zth,nex,ney,nx,ny,nz);

for i_positional_element=1:1:8
    if elements(i_positional_element)~=0
        orientation=9-i_positional_element;
        for ix=1:1:3
            for iy=1:1:3
                for iz=1:1:3
                    
                    
                    solution=solution+w(ix)*w(iy)*w(iz)*(...
                        ...
                        ...
                        terms(elements(i_positional_element),5,ix,iy,iz)...
                        *weights(orientation,type,1,ix,iy,iz)...
                        ...
                        ...
                        +terms(elements(i_positional_element),1,ix,iy,iz)...
                        *weights(orientation,type,2,ix,iy,iz)...
                        +terms(elements(i_positional_element),2,ix,iy,iz)...
                        *weights(orientation,type,3,ix,iy,iz)...
                        +terms(elements(i_positional_element),3,ix,iy,iz)...
                        *weights(orientation,type,4,ix,iy,iz)...
                        ...
                        ...
                        +terms(elements(i_positional_element),4,ix,iy,iz)...
                        *wTerms(orientation,type,ix,iy,iz)...
                        ...
                        ...
                        );
                    
                    
                end
            end
        end
    end
end
solution=solution*dxyz;
end