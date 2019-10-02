function solution=generateWeights()

gp=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
solution=zeros(3,3,16,5);
orientations=...
    [0 0
    0 1
    1 0
    1 1];
types=...
    [0 0
    1 0
    0 1
    1 1];
orders=...
    [0 0
    1 0
    0 1
    2 0
    0 2];

for iorder=1:1:5
    order=orders(iorder,:);
    ilocal=0;
    for iorientation=1:1:4
        for itype=1:1:4
            ilocal=ilocal+1;
            orientation=orientations(iorientation,:);
            type=types(itype,:);
            for i=1:1:3
                solution(i,:,ilocal,iorder)=...
                    basis(gp,orientation(1),type(1),order(1))*...
                    basis(gp(i),orientation(2),type(2),order(2));
            end
        end
    end
end

end