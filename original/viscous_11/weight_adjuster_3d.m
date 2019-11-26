function solution=weight_adjuster_3d(weights,dx,dy,dz)

solution=weights;

types=...
    [0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    1 0 1
    0 1 1
    1 1 1];
orders=...
    [0 0 0
    1 0 0
    0 1 0
    0 0 1
    2 0 0
    0 2 0
    0 0 2];

dhs=[dx dy dz];

for iorientation=1:1:8
    for itype=1:1:8
        for iorder=1:1:7
            for idh=1:1:3
                solution(iorientation,itype,iorder,:,:,:)=...
                    scaling_dh(solution(iorientation,itype,iorder,:,:,:),...
                    dhs(idh),types(itype,idh),orders(iorder,idh));
            end
            
            
%             for igpz=1:1:3
%                 for igpy=1:1:3
%                     for igpx=1:1:3
%                     
%                         
% %                     weights(iorientation,itype,iorder,:,igpy,igpz)=...
% %                         basis(gps,orientations(iorientation,1),types(itype,1),orders(iorder,1))...
% %                         *basis(gps(igpy),orientations(iorientation,2),types(itype,2),orders(iorder,2))...
% %                         *basis(gps(igpz),orientations(iorientation,3),types(itype,3),orders(iorder,3));
%                     
%                     end
%                 end
%             end
        end
    end
end

end

function solution=scaling_dh(vals,dh,type,order)

if type==0
    if order==0
        solution=vals;
    elseif order==1
        solution=vals/dh;
    elseif order==2
        solution=vals/(dh^2);
    end
elseif type==1
    if order==0
        solution=vals*dh;
    elseif order==1
        solution=vals;
    elseif order==2
        solution=vals/dh;
    end
end

end