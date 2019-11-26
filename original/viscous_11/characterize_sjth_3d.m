function sj_properties=characterize_sjth_3d(nnz_,gbfs1,gbfs2,nex,ney,nx,ny,nz,diffT)

%This creates very large matrix. Just passing sj_properties to child
%functions take up too much time! Thus, this method will be disregarded.

if isscalar(diffT)
    sj_properties=zeros(nnz_,8,5);
    parfor sjth=1:1:nnz_
        sj_properties(sjth,:,:)=characterize_sjth_sj_3d(sjth,gbfs1,gbfs2,nex,ney,nx,ny,nz);
    end
else
    sj_properties=zeros(nnz_,8,6);
    parfor sjth=1:1:nnz_
        sj_properties(sjth,:,:)=...
            characterize_sjth_sj_with_gradient_3d(sjth,gbfs1,gbfs2,nex,ney,nx,ny,nz);
    end
end
end

function sjth_properties=characterize_sjth_sj_3d(sjth,gbfs1,gbfs2,nex,ney,nx,ny,nz)
sjth_properties=zeros(8,5);
[elements,lbfs1,lbfs2]=analyze_interaction_3d(...
	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);
for i_positional_element=1:1:8
    if elements(i_positional_element)~=0
        [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
        [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
        sjth_properties(i_positional_element,:)=...
            [elements(i_positional_element),orientation1,orientation2,...
            type1,type2];
    end
end
end

function sjth_properties=characterize_sjth_sj_with_gradient_3d(sjth,gbfs1,gbfs2,nex,ney,nx,ny,nz)
sjth_properties=zeros(8,6);
[elements,lbfs1,lbfs2]=analyze_interaction_3d(...
	gbfs1(sjth),gbfs2(sjth),nex,ney,nx,ny,nz);
for i_positional_element=1:1:8
    if elements(i_positional_element)~=0
        [orientation1,type1]=analyze_lbf_3d(lbfs1(i_positional_element));
        [orientation2,type2]=analyze_lbf_3d(lbfs2(i_positional_element));
        n1xth=get_n1xth_3d(elements(i_positional_element),nex);
        sjth_properties(i_positional_element,:)=...
            [elements(i_positional_element),orientation1,orientation2,...
            type1,type2,n1xth];
    end
end
end