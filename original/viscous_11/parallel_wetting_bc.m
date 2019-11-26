function [sf,sj]=parallel_bc(sj_,iZeros,iOnes,...
    sf_,ll,li,lu,bl,bi,bu,rl,ri,ru,tl,ti,tu,...
    wetting,ig)
sj=sj_;
sf=sf_;

if wetting==0
    for i=iZeros(1,1):iZeros(1,2):iZeros(1,3)
        sj(i)=0.0;
    end
    for i=iZeros(2,1):iZeros(2,2):iZeros(2,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
        sj(i)=0.0;
    end
    for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
        for j=0:1:47
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
        for j=0:1:47
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(7,1):iZeros(7,2):iZeros(7,3)
        sj(i)=0.0;
    end
    for i=iZeros(8,1):iZeros(8,2):iZeros(8,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
        sj(i)=0.0;
    end
    for i=iOnes
        sj(i)=1.0;
    end
    %---Switch on sf-------

    for i=ll:li:lu
        sf(i)=0.0;
        sf(i+2)=0.0;
    end
    for i=bl:bi:bu
        sf(i)=0.0;
        sf(i+1)=0.0;
    end
    for i=rl:ri:ru
        sf(i)=0.0;
        sf(i+2)=0.0;
    end
    for i=tl:ti:tu
        sf(i)=0.0;
        sf(i+1)=0.0;
    end
    
    
elseif wetting==1
    
    
    for i=iZeros(1,1):iZeros(1,2):iZeros(1,3)
        sj(i)=0.0;
    end
    for i=iZeros(2,1):iZeros(2,2):iZeros(2,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
        sj(i)=0.0;
    end
    for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
        sj(i)=0.0;
    end
    for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
        sj(i)=0.0;
    end
    for i=iZeros(7,1):iZeros(7,2):iZeros(7,3)
        sj(i)=0.0;
    end
    for i=iZeros(8,1):iZeros(8,2):iZeros(8,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
        sj(i)=0.0;
    end
    for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
        sj(i)=0.0;
    end
    for i=iZeros(11,1):iZeros(11,2):iZeros(11,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(12,1):iZeros(10,2):iZeros(12,3)
        sj(i)=0.0;
    end
    
    %==================================================
    for i=iZeros(13,1):iZeros(13,2):iZeros(13,3)
        sj(i)=0.0;
    end
    for i=iZeros(14,1):iZeros(14,2):iZeros(14,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(15,1):iZeros(15,2):iZeros(15,3)
        sj(i)=0.0;
    end
    for i=iZeros(16,1):iZeros(16,2):iZeros(16,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(17,1):iZeros(17,2):iZeros(17,3)
        sj(i)=0.0;
    end
    for i=iZeros(18,1):iZeros(18,2):iZeros(18,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    for i=iZeros(19,1):iZeros(19,2):iZeros(19,3)
        sj(i)=0.0;
    end
    for i=iZeros(20,1):iZeros(20,2):iZeros(20,3)
        for j=0:1:23
            sj(i+j)=0.0;
        end
    end
    
    for i=iOnes
        sj(i)=1.0;
    end
    
    [m,~]=size(ig);
    for i=1:1:m
        sj(ig(i,1))=ig(i,2);
    end
    
    for i=ll:li:lu
        sf(i)=0.0;
        sf(i+2)=0.0;
    end
    for i=bl:bi:bu
        sf(i)=0.0;
        sf(i+1)=0.0;
    end
    for i=rl:ri:ru
        sf(i)=0.0;
        sf(i+2)=0.0;
    end
    for i=tl:ti:tu
        sf(i)=0.0;
        sf(i+1)=0.0;
    end
    
%     return
    
%     for i=iZeros(1,1):iZeros(1,2):iZeros(1,3)
%         sj(i)=0.0;
%     end
%     for i=iZeros(2,1):iZeros(2,2):iZeros(2,3)
%         for j=0:1:23
%             sj(i+j)=0.0;
%         end
%     end
% %     for i=iZeros(3,1):iZeros(3,2):iZeros(3,3)
% %         for j=0:1:23
% %             sj(i+j)=0.0;
% %         end
% %     end
%     for i=iZeros(4,1):iZeros(4,2):iZeros(4,3)
%         sj(i)=0.0;
%     end
%     for i=iZeros(5,1):iZeros(5,2):iZeros(5,3)
%         for j=0:1:47
%             sj(i+j)=0.0;
%         end
%     end
% %     for i=iZeros(6,1):iZeros(6,2):iZeros(6,3)
% %         for j=0:1:47
% %             sj(i+j)=0.0;
% %         end
% %     end
%     for i=iZeros(7,1):iZeros(7,2):iZeros(7,3)
%         sj(i)=0.0;
%     end
%     for i=iZeros(8,1):iZeros(8,2):iZeros(8,3)
%         for j=0:1:23
%             sj(i+j)=0.0;
%         end
%     end
% %     for i=iZeros(9,1):iZeros(9,2):iZeros(9,3)
% %         for j=0:1:23
% %             sj(i+j)=0.0;
% %         end
% %     end
%     for i=iZeros(10,1):iZeros(10,2):iZeros(10,3)
%         sj(i)=0.0;
%     end
%     for i=iOnes
%         sj(i)=1.0;
%     end
%     %---Switch on sf-------
% 
%     for i=ll:li:lu
%         sf(i)=0.0;
% %         sf(i+2)=0.0;
%     end
%     for i=bl:bi:bu
%         sf(i)=0.0;
% %         sf(i+1)=0.0;
%     end
%     for i=rl:ri:ru
%         sf(i)=0.0;
% %         sf(i+2)=0.0;
%     end
%     for i=tl:ti:tu
%         sf(i)=0.0;
% %         sf(i+1)=0.0;
%     end
end

end