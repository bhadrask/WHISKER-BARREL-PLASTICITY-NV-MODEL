function[ PO2tissv]=project2vasc(bndry,PO2tis,PO2tis4art,W,loc)
eA=bndry(5); eV=bndry(6);
nE=length(PO2tis);

caps=sqrt(eV-eA-1);

PO2tissue=PO2tis(eA+1:eV-1);
PO2tis2dn=reshape(PO2tissue,caps,caps);

PO2v2d=conv2(PO2tis2dn,W,'same');
PO2v1d=reshape(PO2v2d,caps*caps,1);
PO2v(loc)=PO2v1d;
 po2max=PO2tis4art;
 po2mean=mean(PO2v1d);
 PO2tissv=zeros(nE,1);
 PO2tissv(1:eA)=po2max;
 PO2tissv(eA+1:eV-1)=PO2v;
 PO2tissv(eV:nE)=po2mean/1.5;
   