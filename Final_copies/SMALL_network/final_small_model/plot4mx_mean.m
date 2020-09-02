function plot4mx_mean(x,E,N1,N2,bndry,para,neu_size,nE,L0,S0,Fin0,Fout0,V0,input1,lissom,opts)
k=size(x,1);

%% temporary arrangement
file_vessel='XY.mat'; 
PO2tis=15*ones(nE,1); 
PO2tis4art=30;
W_gauss=gaussian_filter_wt(5,2);
loc_vas=findloc(file_vessel);
kw=3.506*1e-18;
%% 
% if k>20
% nm=1:5:k;
% else
     nm=1:k;
% end
load('XY.mat');
load('vessel_charac.mat'); 
%  load('lissomtrain1.mat'); 
neu_length=neu_size(1)*neu_size(2);
eA=bndry(5); eV=bndry(6);
nE=length(E); nN=N2(end)+1; eA=bndry(5); eV=bndry(6);
caps=eV-eA-1;
d0=0.6*d0;
g=-0.1;

 nA=sum(layers(1:3));
nC=layers(4);
nV=sum(layers(5:7));
S00=x(1,nE+1:2*nE)';
   PO2vess0=exp(0.385*log((1./S00-1).^-1)+3.32-(72.*S00).^-1-(S00.^6)./6);
    PO2tisv0=project2vasc(bndry,PO2tis,PO2tis4art,W_gauss,loc_vas);
        OE0=1*(kw*L0.*d0).*(PO2vess0- PO2tisv0);
beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nA+nC);
CMRO2n0=x(1,3*nE+neu_length+nA+nC+1:end);
lp=1;
ATPref=2.5;
HbOzero=2.3*S0.*V0;
HbTzero=2.3*V0;
Hbzero=2.3*(1-S0).*V0;
meanHbOz=find_mean_ROI(HbOzero(find(lev==4)),4,0,2,6);
meanHbTz=find_mean_ROI(HbTzero(find(lev==4)),4,0,2,6);
meanHbz=find_mean_ROI(Hbzero(find(lev==4)),4,0,2,6);
% pos=[1,2,3,4];
for j=1:length(nm)
   
    g=g+0.1;
    i=nm(j);
V=x(i,1:nE)';
newd=2*sqrt((V*1e9)./(pi*L0));
 R=128*15*(L0/2)*1e3./(pi*newd.^4);

S=x(i,nE+1:2*nE)';
ATP=x(i,3*nE+1:3*nE+neu_length);
% beta=x(i,3*nE+neu_length+1:end);
beta=x(i,3*nE+neu_length+1:3*nE+neu_length+nA+nC);
CMRO2n=x(i,3*nE+neu_length+nA+nC+1:end);
dVdt=zeros(size(V));
PO2tissn1=x(i,2*nE+1:3*nE);
PO2tissn=PO2tissn1(eA+1:eV-1);
% [Pn, Pe, Fin, Fout]= vsc_func6(t,R,V,E,N1,N2,bndry,para,beta); %Finding Flow rate

[Pnv,Pe]=vsc_func2_old(E,N1,N2,bndry,R,V,para,beta);
[PnA,Pe]=vsc_func3(E,N1,N2,bndry,R,Pe,dVdt,para);
Pn=[PnA;Pnv];
[dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe);
Fin1=zeros(2,size(Fin,1));
Fout1=zeros(2,size(Fin,1));
Fin1(1,:)=Fin0'; Fout1(1,:)=Fout0';
Fin1(2,:)=Fin';Fout1(2,:)=Fout';
F=(Fin1+Fout1)/2;
V1=zeros(2,size(Fin0,1));
S1=zeros(2,size(Fin0,1));
V1(1,:)=V0';
V1(2,:)=V';
S1(1,:)=S0';
S1(2,:)=S';

 HbT=2.3*V1;
 HbO=2.3.*S1.*V1;
 Hb=2.3*(1-S1).*V1;
 BETA(1,:)=beta0';
 BETA(2,:)=beta';
 if g>1 && g<=2
   
    input=input1;
 else
    input=zeros(size(input1,1),size(input1,2));
 end

[neu,lissom]=inst_response_updat_depending_vsclr_fdk(lissom,opts,input,ATP,0);
% if g>=4.8

% end
meanHbO=find_mean_ROI(HbO(2,find(lev==4)),4,0,2,6);
meanHbT=find_mean_ROI(HbT(2,find(lev==4)),4,0,2,6);
meanHb=find_mean_ROI(Hb(2,find(lev==4)),4,0,2,6);

    PO2vess=exp(0.385*log((1./S-1).^-1)+3.32-(72.*S).^-1-(S.^6)./6);
    PO2tisv=project2vasc(bndry,PO2tissn1,PO2tis4art,W_gauss,loc_vas);
        OEf=(((kw*L0.*newd.*(PO2vess- PO2tisv))-OE0)./OE0);
        OE(lp)=1+mean(OEf(find(lev==4)));
 meancmro2(lp)=1+mean((CMRO2n-CMRO2n0)./CMRO2n0);
        meanatp(lp)=mean(ATP);
        meanF(lp)=1+mean((F(2,eA+1:eV-1)-F(1,eA+1:eV-1))./F(1,eA+1:eV-1));
        meanpo2(lp)=1+mean((PO2tissn-15*ones(size(PO2tissn)))./(15*ones(size(PO2tissn))));
       meansat(lp)=1+mean((S1(2,eA+1:eV-1)-S1(1,find(lev==4)))./S1(1,find(lev==4)));

       meanbeta(lp)=mean(beta(eA+1:eV-1));
       MHbT(lp)=1+(meanHbO-meanHbOz)/meanHbOz;
       MHbO(lp)=1+(meanHbT-meanHbTz)/meanHbTz;
       MHb(lp)=1+(meanHb-meanHbz)/meanHbz;
       lp=lp+1;


end
T=0:0.1:5;
figure();plot(T,meancmro2);title(['cmro2  ']);%,num2str((max(meancmro2)-min(meancmro2))*100/min(meancmro2))]); hold on; plot(T,repmat(mean(CMRO2nref),size(meancmro2)))
figure();plot(T,meanatp);title(['atp  ']);%,num2str((max(meanatp)-min(meanatp))*100/min(meanatp))]);
figure();plot(T,meanF);title(['Flow  ',num2str((max(meanF)-min(meanF))*100/min(meanF))]);
figure();plot(T,meanpo2);title(['PO2 ',num2str((max(meanpo2)-min(meanpo2))*100/min(meanpo2))]);
% figure();plot(T,meanoe);title(['OEn ',num2str((max(meanoe)-min(meanoe))*100/min(meanoe))]);
% figure();plot(T,meanpo2vess);title(['PO2_vessel',num2str((max( meanpo2vess)-min( meanpo2vess))*100/min( meanpo2vess))]);
figure();plot(T,meansat);title(['S',num2str((max( meansat)-min( meansat))*100/min( meansat))]);
% figure();plot(T,meanN);title(['N']);
figure();plot(T,meanbeta);title(['BETA']);
figure();plot(T,MHbT);title(['HbT',num2str((max(MHbT)-min(MHbT))*100/min(MHbT))]);
figure();plot(T,MHbO);title(['HbO',num2str((max(MHbO)-min(MHbO))*100/min(MHbO))]);
figure();plot(T,MHb);title(['Hb',num2str((max(MHb)-min(MHb))*100/min(MHb))]);
figure();plot(T,OE);title('OE');
figure();plot(T,meancmro2);hold on;plot(T,MHbT,'r') ;hold on;plot(T,MHbO,'g');