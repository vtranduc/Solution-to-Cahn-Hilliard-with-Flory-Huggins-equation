function weights=set_up_rec_mesh_3D(gp)

weights=zeros(64,7,27);

k=0;
for w1=1:1:3
    for w2=1:1:3
        for w3=1:1:3
            k=k+1;
            [phi,phix,phiy,phiz,phixx,phiyy,phizz]=get_basis(gp(w1),gp(w2),gp(w3));
            weights(:,:,k)=[phi',phix',phiy',phiz',phixx',phiyy',phizz'];
        end
    end
end

end