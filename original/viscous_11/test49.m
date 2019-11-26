function test49(X,Y,Z,diffT,chi,nex,nx,ny,nexney,ne,T,diff)

slope=T(2)-T(1);

gps=[(-sqrt(3/5)+1)/2 0.5 (sqrt(3/5)+1)/2];
ne
eTest=10
plot3(X,Y,Z,'y.')
grid on
xlabel('x');ylabel('y');zlabel('z')
for e=1:1:ne
    [eX,eY,eZ]=get_eXYZ_3d(e,nex,nx,ny,nexney,X,Y,Z);
    
    if e==eTest
        hold on
        plot3(eX,eY,eZ,'ro')
    end
    
    for i=1:1:8
        dist=dist_3d(eX(i),eY(i),eZ(i));
        if dist>1
            dist
            error('hello')
        end
        
        for igamma=1:1:3
            for ibeta=1:1:3
                for ialpha=1:1:3
                    x=map_to_global_1_compt_3d(gps(ialpha),gps(ibeta),gps(igamma),eX,[0 0 0]);
                    y=map_to_global_1_compt_3d(gps(ialpha),gps(ibeta),gps(igamma),eY,[0 0 0]);
                    z=map_to_global_1_compt_3d(gps(ialpha),gps(ibeta),gps(igamma),eZ,[0 0 0]);
                    
                    if e==eTest
                        hold on
                        plot3(x,y,z,'g.')
                        hold off
                    end
                    
                    dist=sqrt(x^2+y^2+z^2);
                    
                    err=abs((dist*slope+T(1))-diffT(e,ialpha,ibeta,igamma)/diff);
                    
                    if err>10^-16
                        dist_3d(eX(1),eY(1),eZ(1))
                        [eX' eY' eZ']
                        plot3(X,Y,Z,'r.',x,y,z,'bo')
                        grid on
                        e
                        [x y z]
                        dist
                        slope
                        diffT(e,ialpha,ibeta,igamma)/diff
                        err
                        error('bad')
                    end
                    
                end
            end
        end
        
    end
end



end

function dist=dist_3d(x,y,z)
dist=sqrt(x^2+y^2+z^2);
end