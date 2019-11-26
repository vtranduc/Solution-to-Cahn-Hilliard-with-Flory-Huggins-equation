function gps_transformed_mx=get_transformed_weight_points_3d(xyz_mx,bases)

%==== THIS SHOULD BE PASSED AS ARGUMENT
% gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
% w_1comp=zeros(2,3);
% w_1comp(2,:)=gps;
% for i=1:1:3
%     w_1comp(1,i)=1-gps(i);
% end
% bases=zeros(3,3,3,8);
% for gammath=1:1:3
%     for betath=1:1:3
%         for alphath=1:1:3
%             bases(alphath,betath,gammath,1)=w_1comp(1,alphath)*w_1comp(1,betath)*w_1comp(1,gammath);
%             bases(alphath,betath,gammath,2)=w_1comp(2,alphath)*w_1comp(1,betath)*w_1comp(1,gammath);
%             bases(alphath,betath,gammath,3)=w_1comp(1,alphath)*w_1comp(2,betath)*w_1comp(1,gammath);
%             bases(alphath,betath,gammath,4)=w_1comp(2,alphath)*w_1comp(2,betath)*w_1comp(1,gammath);
%             bases(alphath,betath,gammath,5)=w_1comp(1,alphath)*w_1comp(1,betath)*w_1comp(2,gammath);
%             bases(alphath,betath,gammath,6)=w_1comp(2,alphath)*w_1comp(1,betath)*w_1comp(2,gammath);
%             bases(alphath,betath,gammath,7)=w_1comp(1,alphath)*w_1comp(2,betath)*w_1comp(2,gammath);
%             bases(alphath,betath,gammath,8)=w_1comp(2,alphath)*w_1comp(2,betath)*w_1comp(2,gammath);
%         end
%     end
% end
%==================

gps_transformed_mx=zeros(3,3,3,3);

for gammath=1:1:3
    for betath=1:1:3
        for alphath=1:1:3
            for compth=1:1:3
                pos_new=0;
                for edgeth=1:1:8
                    pos_new=pos_new+xyz_mx(edgeth,compth)*bases(alphath,betath,gammath,edgeth);
                end
                gps_transformed_mx(alphath,betath,gammath,compth)=pos_new;
            end
        end
    end
end

end