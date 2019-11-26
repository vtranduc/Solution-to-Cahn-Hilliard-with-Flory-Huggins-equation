function solution=tfunct3D(vals,alpha)

%Vals must be in vector form, containing 3 elements

orientation=[0 0 0
    0 1 0
    1 0 0
    1 1 0
    0 0 1
    0 1 1
    1 0 1
    1 1 1];

type=[0 0 0
    1 0 0
    0 1 0
    0 0 1
    1 1 0
    1 0 1
    0 1 1
    1 1 1];

order=[0 0 0
    1 0 0
    0 1 0
    0 0 1
    2 0 0
    0 2 0
    0 0 2];

solution=zeros(64,7);

for iorder=1:1:7
    k=0;
    for inode=1:1:8
        for itype=1:1:8
            k=k+1;
            phi_=1;
            for icomp=1:1:3
                phi_=phi_*basis_(vals(icomp),orientation(inode,icomp),...
                    type(itype,icomp),order(iorder,icomp),alpha(icomp));
            end
%             if phi_==0
%                 vals
%                 orientation(inode,icomp)
%                 type(itype,icomp)
%                 order(iorder,icomp)
%                 for icomp=1:1:3
%                 	basis(vals(icomp),orientation(inode,icomp),type(itype,icomp),order(iorder,icomp))
%                 end
%                 error('ZERO DETECTED')
%             end
            solution(k,iorder)=phi_;
        end
    end
end

end