function sj=test66(sjr,neight,n_nSurf,neightTotal,inverse_rotation)

sj=sjr;
first_index=neight-7;
for i=1:1:n_nSurf(6)
    first_index=first_index+8;
    
    for typer=1:1:3
        
        
        for col=neight:1:neightTotal
            [~,type]=analyze_gbs_3d(col);
            if type==2
                sj(first_index+typer,col)=...
                    inverse_rotation(i,1,1)*sjr(first_index+typer,col)...
                    +inverse_rotation(i,2,1)*sjr(first_index+typer,col+1)...
                    +inverse_rotation(i,3,1)*sjr(first_index+typer,col+2);
            end
            if type==3
                sj(first_index+typer,col)=...
                    inverse_rotation(i,1,2)*sjr(first_index+typer,col-1)...
                    +inverse_rotation(i,2,2)*sjr(first_index+typer,col)...
                    +inverse_rotation(i,3,2)*sjr(first_index+typer,col+1);
            end
            if type==4
                sj(first_index+typer,col)=...
                    inverse_rotation(i,1,3)*sjr(first_index+typer,col-2)...
                    +inverse_rotation(i,2,3)*sjr(first_index+typer,col-1)...
                    +inverse_rotation(i,3,3)*sjr(first_index+typer,col);
            end
        end
    
    end
    
end


% sj=sjr;
% first_index=neight-7;
% for i=1:1:n_nSurf(6)
%     first_index=first_index+8;
%     for ijk=1:1:neightTotal
% %         if sjr(first_index+1,ijk)~=0 %-----------
%         if 1==1
%             for j=1:1:3
%                 sj(first_index+j,ijk)=...
%                     inverse_rotation(i,j,1)*sjr(first_index+1,ijk)...
%                     +inverse_rotation(i,j,2)*sjr(first_index+2,ijk)...
%                     +inverse_rotation(i,j,3)*sjr(first_index+3,ijk);
%             end
%         end
%     end
% end

end