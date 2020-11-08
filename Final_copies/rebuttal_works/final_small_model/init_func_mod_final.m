function[ OE,OEn,CMRO2n]=init_func_mod_final(S,bndry,neu_size,PO2tissv,L0,d,W,filename,file_vessel_chara,OEn0)

eA=bndry(5); eV=bndry(6);
 kw=3.506*1e-18;
 OE=zeros(1,length(S));
 load(filename);
 load(file_vessel_chara,'layers');
 cap_layer=(numel(layers)+1)/2;

caps=sqrt(eV-eA-1);

Gx=X(find(lev==cap_layer));
Gy=Y(find(lev==cap_layer));

 
 for i=1:length(S)
    PO2=exp(0.385*log((1/S(i)-1)^-1)+3.32-(72*S(i))^-1-(S(i)^6)/6);
    po2(i)=PO2;
   
    OE(i)=1*(kw*L0(i)*d(i))*(PO2-PO2tissv(i));
   
   
 end
 OE(OE<0)=0;
 oe2d=zeros(caps,caps);
OE2=OE(eA+1:eV-1);
  count=0;
for i=1:length(OE2)
    oe2d(Gx(i),Gy(i))=OE2(i);
    
end
  OEn=zeros(neu_size(1)*neu_size(2),1);
   CMRO2n=zeros(neu_size(1)*neu_size(2),1);
   oe2dn=conv2(oe2d,W,'same');
  
 for cj=1:neu_size(1)
     for ci=1:neu_size(2)

count=count+1;


OEn(count)=oe2dn(ci,cj);
CMRO2n(count)= OEn0(count);
     end
 end
 