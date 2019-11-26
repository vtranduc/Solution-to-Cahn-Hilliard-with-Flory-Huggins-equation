function sol=map_to_global_1_compt_2d(alpha,beta,eCompt,orders)
sol=eCompt(1)*isoBasis(alpha,0,orders(1))*isoBasis(beta,0,orders(2))...
    +eCompt(2)*isoBasis(alpha,1,orders(1))*isoBasis(beta,0,orders(2))...
    +eCompt(3)*isoBasis(alpha,0,orders(1))*isoBasis(beta,1,orders(2))...
    +eCompt(4)*isoBasis(alpha,1,orders(1))*isoBasis(beta,1,orders(2));
end