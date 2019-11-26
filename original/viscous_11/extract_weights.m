function [phi,phix,phiy,phiz,phixx,phiyy,phizz]=extract_weights(iweight,weights)

phi=weights(:,1,iweight);
phix=weights(:,2,iweight);
phiy=weights(:,3,iweight);
phiz=weights(:,4,iweight);
phixx=weights(:,5,iweight);
phiyy=weights(:,6,iweight);
phizz=weights(:,7,iweight);

% phi=phi';
% phix=phix';
% phiy=phiy';
% phiz=phiz';
% phixx=phixx';
% phiyy=phiyy';
% phizz=phizz';

end