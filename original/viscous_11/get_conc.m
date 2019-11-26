function solution=get_conc(ne,ney,ny,c,weights)

solution=zeros(3,3,5,ne);

parfor element=1:1:ne
    [inx,iny]=inxiny_elemental(element,ney);
    gbfs=elemental_gbf(inx,iny,ny);
    cElemental=get_cElemental(gbfs,c);
    
    
    conc_=zeros(3,3,5);
    for order_type=1:1:5
        conc_(:,:,order_type)=conc(cElemental,weights,order_type);
    end
    
    solution(:,:,:,element)=conc_;
    
end

end