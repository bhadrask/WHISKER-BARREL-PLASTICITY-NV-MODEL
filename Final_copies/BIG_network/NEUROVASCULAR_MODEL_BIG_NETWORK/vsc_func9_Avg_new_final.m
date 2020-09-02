function [dATPdt,Jatp]=vsc_func9_Avg_new_final(CMRO2n,ATP,neu,lissom,CMRO2nref,ATPref)

neu1=reshape(neu,lissom.layers{2}.dim(1)* lissom.layers{2}.dim(2),1);

factr=1;
term1=(ATPref.*factr-ATP);
term2=((CMRO2n-CMRO2nref)./CMRO2nref);
term2(isnan(term2))=0;
term2(isinf(term2))=0;
term3=term2-0.5*neu1; %"Nat"

dATPdt=(term1+term3)/4; %tatp

Jatp=term2;
end