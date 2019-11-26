function solution=dim_reducer1_3d(vals)
solution=zeros(3,3,3);
for iz=1:1:3
    for iy=1:1:3
        for ix=1:1:3
            solution(ix,iy,iz)=vals(1,ix,iy,iz);
        end
    end
end
end