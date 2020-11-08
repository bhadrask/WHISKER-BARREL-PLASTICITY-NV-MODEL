function [dCMRO2ndt,OEn]=vsc_func7_Avg_bsk_mod(OE1,bndry,neu,CMRO2nref,W,filename,CMRO2n,PO2tis,PO2tismax,layers)
eA=bndry(5); eV=bndry(6);
cap_layer=(numel(layers)+1)/2;
neu_size=[size(neu,1),size(neu,2)];
OE2=OE1(eA+1:eV-1);
load(filename);
G=[X(find(lev==cap_layer)),Y(find((lev==cap_layer)))];
eA=bndry(5); eV=bndry(6);
caps=sqrt(eV-eA-1);
oe2d=zeros(caps,caps);
for i=1:length(OE2)
    oe2d(G(i,1),G(i,2))=OE2(i);
end
% OEn=zeros(neu_size(1)*neu_size(2),1);
% dCMRO2ndt=zeros(neu_size(1)*neu_size(2),1);
% count=0;
oe2dn=conv2(oe2d,W,'same');
OEn=reshape(oe2dn,neu_size(1)*neu_size(2),1);
PO2tissue=PO2tis(eA+1:eV-1);
PO2tissue_max=PO2tismax(eA+1:eV-1);
 factr=PO2tissue./PO2tissue_max;
  factr(factr>1)=1;
  dCMRO2ndt=(0.7*OEn.*(neu)+(CMRO2nref.*factr-CMRO2n))/2;

% for cj=1:size(neu,1)
%     for ci=1:size(neu,2)
%         count=count+1;
%         OEn(count)=oe2dn(ci,cj);
%         dCMRO2ndt(count)=(0.7*OEn(count)*(neu(ci,cj))+(CMRO2nref(count)*factr(count)-CMRO2n(count)))/2;
% %                 dCMRO2ndt(count)=(0.001*OEn(count)*neu(ci,cj)+(CMRO2nref(count)-CMRO2n(count)))/0.1;
% %                 older working one
%     end
% end
end