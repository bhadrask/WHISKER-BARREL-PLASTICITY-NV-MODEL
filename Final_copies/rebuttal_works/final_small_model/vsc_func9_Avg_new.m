function [dATPdt,Jatp]=vsc_func9_Avg_new(CMRO2n,ATP,neu1,lissom,CMRO2nref,ATPref,PO2tis,PO2tismax,bndry)
tau=1e-2;
% neu1=reshape(neu,lissom.layers{2}.dim(1)* lissom.layers{2}.dim(2),1);
% locs=find(OEn>OE0);
% OEtr=OEn; OEtr(locs)=OE0(locs);
%      keyboard;
%  eA=bndry(5); eV=bndry(6);

% PO2tissue=PO2tis(eA+1:eV-1);
% PO2tissue_max=PO2tismax(eA+1:eV-1);
%  factr=PO2tissue./PO2tissue_max;
%   factr(factr>1)=1;
%  factr(factr<0.9999)=0.8;
factr=1;
term1=(ATPref.*factr-ATP);
term2=((CMRO2n-CMRO2nref)./CMRO2nref);
term2(isnan(term2))=0;
term2(isinf(term2))=0;
term3=term2-0.5*neu1; %"Nat"

dATPdt=(term1+term3)/4; %tatp

Jatp=term2;
end