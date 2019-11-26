function sj=specific_row_reset_rec_3D(sj,row,ny,nz,nxny)

type=mod(row,8);
if type==0
    type=8;
end

node=(row-type)/8+1;

interactions=interactions_rec_3D(node,ny,nz,nxny);
for k=1:1:length(interactions)
	pre_index=interactions(k)*8-8;
%     interactions(k)
%     fprintf('is the interaction')
%     node
	for l=1:1:8
%         row
%         pre_index+l
%         if sj(row,pre_index+l)==0
%             error('just raise')
%         end
        sj(row,pre_index+l)=0.0;
	end
end


end