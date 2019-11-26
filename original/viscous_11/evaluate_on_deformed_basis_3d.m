function evaluate_on_deformed_basis_3d(points,xyz_mx)

% clear
% clc



[rt1,rt2,rt3,len]=rotate_arm_on_x()



end

function eval_rotated_basis(pt,rt1,rt2,rt3,len,orientation,type,order)

term=pt(1)*rt1-pt(2)*rt2+pt(3)*rt3;
basis(term,orientation,type,order);

end

function [rt1,rt2,rt3,len]=rotate_arm_on_x(xyz_mx,elemental_basis)

%-------------- TEST
start=zeros(1,3); %NEXT, WE SHOULD CONSIDER THE NONZERO CASE!!!
tip=[-7,13,-sqrt(7)];
%--------------First rotation

[unit_vector,len]=get_unit_vector(start,tip);

theta1=-arctansp(unit_vector(1),unit_vector(2));
adj_len=sqrt(unit_vector(1)^2+unit_vector(2)^2);
theta2=-arctansp(adj_len,unit_vector(3));

sin1=sin(theta1);
sin2=sin(theta2);
cos1=cos(theta1);
cos2=cos(theta2);

rt1=cos1;
rt2=sin1*cos2;
rt3=sin1*sin2;

end

function [unit_vector,len]=get_unit_vector(start,tip)
% It points from "start" to "tip"
pointer=tip-start;
len=sqrt(sum(pointer.^2));
unit_vector=pointer/len;
end


% 
% function len=get_length(start,tip)
% len=sqrt((tip-start).^2);
% end