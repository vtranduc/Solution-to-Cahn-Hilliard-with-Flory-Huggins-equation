function weights_plus=generateWeights_plus_3d()

weights_plus=zeros(8,8,10,3,3,3);

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];

orientations=...
    [0 0 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1
    1 0 1
    0 1 1
    1 1 1];
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
    0 0 2
    1 1 0
    1 0 1
    0 1 1];

for iorientation=1:1:8
    for itype=1:1:8
        for iorder=1:1:10
            for igpz=1:1:3
                for igpy=1:1:3
                    weights_plus(iorientation,itype,iorder,:,igpy,igpz)=...
                        basis(gps,orientations(iorientation,1),types(itype,1),orders(iorder,1))...
                        *basis(gps(igpy),orientations(iorientation,2),types(itype,2),orders(iorder,2))...
                        *basis(gps(igpz),orientations(iorientation,3),types(itype,3),orders(iorder,3));
                end
            end
        end
    end
end

end