function solution=elemental_gbf(inx_ll,iny_ll,ny)
solution=zeros(1,16);
adder=ny*4;
pre_gbf=adder*(inx_ll-1)+(iny_ll-1)*4;
for i=1:1:8
    solution(i)=pre_gbf+i;
    solution(i+8)=solution(i)+adder;
end
end