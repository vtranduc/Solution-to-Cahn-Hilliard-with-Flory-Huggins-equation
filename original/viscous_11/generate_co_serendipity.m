function co=generate_co_serendipity(ci,include_fluc,ci_fluc,n,nfour,nx,ny,nz)
co=zeros(1,nfour);
if length(ci)==1
    adder=ci-ci_fluc*include_fluc;
    for i=1:1:n
        co(4*(i-1)+1)=adder+2*ci_fluc*include_fluc*rand(1,1);
    end
    deviation=sum(co)/n-ci;
    for i=1:1:n
        co(4*(i-1)+1)=co(4*(i-1)+1)-deviation;
    end
end
% elseif length(ci)==2
%     node=0;
%     adder=linspace(ci(1),ci(2),nx)-ci_fluc*include_fluc*ones(1,nx);
%     for zth=1:1:nz
%         for yth=1:1:ny
%             for xth=1:1:nx
%                 node=node+1;
%                 co(8*(node-1)+1)=adder(xth)+2*ci_fluc*include_fluc*rand(1,1);
%             end
%         end
%     end
%     ci_ave=mean(ci);
%     deviation=sum(co)/n-ci_ave;
%     for i=1:1:n
%         co(8*(i-1)+1)=co(8*(i-1)+1)-deviation;
%     end
% end
end