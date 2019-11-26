function test61

clear
clc

% xyz=[1,-1,-1]
% 
% 
% 
% 
% spherical_coord=cartesian_to_spherical(xyz);
% % 
% % % spherical_coord(2)=spherical_coord(2)*180/pi;
% % % spherical_coord(3)=spherical_coord(3)*180/pi;
% % 
% % 
% % spherical_coord
% % 
% % sqrt(3)
% 
% xyz=spherical_to_cartesian(spherical_coord);

for x=-1:0.1:1
    for y=-1:0.1:1
        for z=-1:0.1:1
            xyz2=[x,y,z];
            spherical_coord=cartesian_to_spherical(xyz2);
            xyz=spherical_to_cartesian(spherical_coord);
            
            
            diff=sum(abs(xyz-xyz2));
            if diff>10^-15
                xyz
                xyz2
                spherical_coord*180/pi
                diff
                error('bad')
            end
            
        end
    end
end

end