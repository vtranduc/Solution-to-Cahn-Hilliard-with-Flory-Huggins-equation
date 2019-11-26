function sol=compute_elemental_weights_serendipity(gps,phis,eX,eY,eZ)
sol=zeros(8,4,7,3,3,3);

ders=zeros(3,7);

derivative_local=zeros(9,1);

for ix=1:1:3
    for iy=1:1:3
        for iz=1:1:3

            ders(1,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 0 0]);
            ders(1,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[0 1 0]);
            ders(1,3)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[0 0 1]);
            ders(1,4)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 1 0]);
            ders(1,5)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 0 1]);
            ders(1,6)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[0 1 1]);
            ders(1,7)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eX,[1 1 1]);

            ders(2,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[1 0 0]);
            ders(2,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 1 0]);
            ders(2,3)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 0 1]);
            ders(2,4)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[1 1 0]);
            ders(2,5)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[1 0 1]);
            ders(2,6)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[0 1 1]);
            ders(2,7)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eY,[1 1 1]);

            ders(3,1)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[1 0 0]);
            ders(3,2)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[0 1 0]);
            ders(3,3)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[0 0 1]);
            ders(3,4)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[1 1 0]);
            ders(3,5)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[1 0 1]);
            ders(3,6)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[0 1 1]);
            ders(3,7)=map_to_global_1_compt_3d(gps(ix),gps(iy),gps(iz),eZ,[1 1 1]);
            
            J=get_isoparametric_mapping_Jacobian(gps(ix),gps(iy),gps(iz),eX,eY,eZ);
            
            for orientation=1:1:8
                
                sol(orientation,1,1,ix,iy,iz)=...
                    phis(orientation,1,1,ix,iy,iz);
                sol(orientation,2,1,ix,iy,iz)=...
                    phis(orientation,2,1,ix,iy,iz)*ders(1,1)...
                    +phis(orientation,3,1,ix,iy,iz)*ders(1,2)...
                    +phis(orientation,4,1,ix,iy,iz)*ders(1,3);
                sol(orientation,3,1,ix,iy,iz)=...
                    phis(orientation,2,1,ix,iy,iz)*ders(2,1)...
                    +phis(orientation,3,1,ix,iy,iz)*ders(2,2)...
                    +phis(orientation,4,1,ix,iy,iz)*ders(2,3);
                sol(orientation,4,1,ix,iy,iz)=...
                    phis(orientation,2,1,ix,iy,iz)*ders(3,1)...
                    +phis(orientation,3,1,ix,iy,iz)*ders(3,2)...
                    +phis(orientation,4,1,ix,iy,iz)*ders(3,3);
                
                for order=2:1:10
                    derivative_local(order-1)=...
                        phis(orientation,1,order,ix,iy,iz);
                end
                x=J\derivative_local;
                for order=2:1:7
                    sol(orientation,1,order,ix,iy,iz)=x(order-1);
                end
                
                for type=1:1:3
                
                    derivative_local(1)=...
                        phis(orientation,2,2,ix,iy,iz)*ders(type,1)+phis(orientation,3,2,ix,iy,iz)*ders(type,2)...
                        +phis(orientation,3,1,ix,iy,iz)*ders(type,4)+phis(orientation,4,2,ix,iy,iz)*ders(type,3)...
                        +phis(orientation,4,1,ix,iy,iz)*ders(type,5);
                    derivative_local(2)=...
                        phis(orientation,2,3,ix,iy,iz)*ders(type,1)+phis(orientation,2,1,ix,iy,iz)*ders(type,4)...
                        +phis(orientation,3,3,ix,iy,iz)*ders(type,2)+phis(orientation,4,3,ix,iy,iz)*ders(type,3)...
                        +phis(orientation,4,1,ix,iy,iz)*ders(type,6);
                    derivative_local(3)=...
                        phis(orientation,2,4,ix,iy,iz)*ders(type,1)+phis(orientation,2,1,ix,iy,iz)*ders(type,5)...
                        +phis(orientation,3,4,ix,iy,iz)*ders(type,2)+phis(orientation,3,1,ix,iy,iz)*ders(type,6)...
                        +phis(orientation,4,4,ix,iy,iz)*ders(type,3);
                    derivative_local(4)=...
                        phis(orientation,2,5,ix,iy,iz)*ders(type,1)+phis(orientation,3,5,ix,iy,iz)*ders(type,2)...
                        +phis(orientation,3,2,ix,iy,iz)*ders(type,4)+phis(orientation,3,2,ix,iy,iz)*ders(type,4)...
                        +phis(orientation,4,5,ix,iy,iz)*ders(type,3)+phis(orientation,4,2,ix,iy,iz)*ders(type,5)...
                        +phis(orientation,4,2,ix,iy,iz)*ders(type,5);
                    derivative_local(5)=...
                        phis(orientation,2,6,ix,iy,iz)*ders(type,1)+phis(orientation,2,3,ix,iy,iz)*ders(type,4)...
                        +phis(orientation,2,3,ix,iy,iz)*ders(type,4)+phis(orientation,3,6,ix,iy,iz)*ders(type,2)...
                        +phis(orientation,4,6,ix,iy,iz)*ders(type,3)+phis(orientation,4,3,ix,iy,iz)*ders(type,6)...
                        +phis(orientation,4,3,ix,iy,iz)*ders(type,6);
                    derivative_local(6)=...
                        phis(orientation,2,7,ix,iy,iz)*ders(type,1)+phis(orientation,2,4,ix,iy,iz)*ders(type,5)...
                        +phis(orientation,2,4,ix,iy,iz)*ders(type,5)+phis(orientation,3,7,ix,iy,iz)*ders(type,2)...
                        +phis(orientation,3,4,ix,iy,iz)*ders(type,6)+phis(orientation,3,4,ix,iy,iz)*ders(type,6)...
                        +phis(orientation,4,7,ix,iy,iz)*ders(type,3);
                    derivative_local(7)=...
                        phis(orientation,2,8,ix,iy,iz)*ders(type,1)+phis(orientation,2,2,ix,iy,iz)*ders(type,4)...
                        +phis(orientation,3,8,ix,iy,iz)*ders(type,2)+phis(orientation,3,3,ix,iy,iz)*ders(type,4)...
                        +phis(orientation,4,8,ix,iy,iz)*ders(type,3)+phis(orientation,4,2,ix,iy,iz)*ders(type,6)...
                        +phis(orientation,4,3,ix,iy,iz)*ders(type,5)+phis(orientation,4,1,ix,iy,iz)*ders(type,7);
                    derivative_local(8)=...
                        phis(orientation,2,9,ix,iy,iz)*ders(type,1)+phis(orientation,2,2,ix,iy,iz)*ders(type,5)...
                        +phis(orientation,3,9,ix,iy,iz)*ders(type,2)+phis(orientation,3,2,ix,iy,iz)*ders(type,6)...
                        +phis(orientation,3,4,ix,iy,iz)*ders(type,4)+phis(orientation,3,1,ix,iy,iz)*ders(type,7)...
                        +phis(orientation,4,9,ix,iy,iz)*ders(type,3)+phis(orientation,4,4,ix,iy,iz)*ders(type,5);
                    derivative_local(9)=...
                        phis(orientation,2,10,ix,iy,iz)*ders(type,1)+phis(orientation,2,3,ix,iy,iz)*ders(type,5)...
                        +phis(orientation,2,4,ix,iy,iz)*ders(type,4)+phis(orientation,2,1,ix,iy,iz)*ders(type,7)...
                        +phis(orientation,3,10,ix,iy,iz)*ders(type,2)+phis(orientation,3,3,ix,iy,iz)*ders(type,6)...
                        +phis(orientation,4,10,ix,iy,iz)*ders(type,3)+phis(orientation,4,4,ix,iy,iz)*ders(type,6);
                    x=J\derivative_local;
                                   
                    for order=1:1:6
                        sol(orientation,type+1,order+1,ix,iy,iz)=x(order);
                    end
                    
                end
            end
        end
    end
end




end