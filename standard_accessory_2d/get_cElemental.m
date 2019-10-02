function solution=get_cElemental(gbfs,c)
solution=zeros(1,16);
for i=1:1:16
    solution(i)=c(gbfs(i));
end
end