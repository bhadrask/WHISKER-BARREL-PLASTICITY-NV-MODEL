 close all;
clc;
clear all;

% filename = 'K1_2cap_2art_K2_0.05_0.5s_cmro2delay_new_pv_reltn_1.7_exponent_0.5s_tau_2Cap_2_art_betadelay.mat';%sprintf('r%d.mat',inps(num)) ;
 filename = 'TESTv1_lesion_sat0.4_atpthr20.mat';
 load(filename);k=size(x,1);
load('XY.mat');
load('vessel_charac.mat'); 
% input1=IN.input;
% if k>20
% nm=1:5:k;
% else
     nm=1:k;
% end

neu_length=neu_size(1)*neu_size(2);
eA=bndry(5); eV=bndry(6);

g=-0.1;


nA=sum(layers(1:3));
nC=layers(4);
nV=sum(layers(5:7));




beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n0=x(1,3*nE+neu_length+nE+1:end);

caps=sqrt(nC);

area_select1=zeros(caps,caps);
area_select1(2:6,1:4)=1;
lp=1;
ATPref=2.5;
loopnr=0;
for j=1:length(nm)
   
  
    i=nm(j);
    loopnr=loopnr+1;
V=x(i,1:nE)';
newd=2*sqrt((V*1e9)./(pi*L0));
 R=128*15*(L0/2)*1e3./(pi*newd.^4);

V=x(i,1:nE);
S=x(i,nE+1:2*nE);
PO2tis=x(i,2*nE+1:3*nE);
ATP=x(i,3*nE+1:3*nE+neu_length);
% ATP_prev=z(3*nE+1:3*nE+neu_length);
beta=x(i,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n=x(i,3*nE+neu_length+nE+1:end);
dVdt=zeros(size(V));
PO2tissn1=x(i,2*nE+1:3*nE);
PO2tissn=PO2tissn1(eA+1:eV-1);
% [Pn, Pe, Fin, Fout]= vsc_func6(t,R,V,E,N1,N2,bndry,para,beta); %Finding Flow rate

[Pnv,Pe,dv]=vsc_func2_new(E,N1,N2,bndry,R,V,para,beta);
[PnA,Pe]=vsc_func3_new(E,N1,N2,bndry,R,Pe,dVdt,para);
Pn=[PnA;Pnv];
[dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe);
if j==1
    Fin0=Fin;
    Fout0=Fout;
end
Fin1=zeros(2,size(Fin,1));
Fout1=zeros(2,size(Fin,1));
Fin1(1,:)=Fin0'; Fout1(1,:)=Fout0';
Fin1(2,:)=Fin';Fout1(2,:)=Fout';
F=(Fin1+Fout1)/2;
V1=zeros(2,size(Fin0,1));
S1=zeros(2,size(Fin0,1));
CMR=zeros(2,length(CMRO2n0));
CMR(1,:)=CMRO2n0';
CMR(2,:)=CMRO2n';
V1(1,:)=V0';
V1(2,:)=V';
S1(1,:)=S0';
S1(2,:)=S';
 HbT=2.3*V1;
 HbO=2.3.*S1.*V1;
 Hb=2.3*(1-S1).*V1;
 BETA(1,:)=beta0';
 BETA(2,:)=beta';

ATP2d=reshape(ATP,8,8);
CMRO2_2d=reshape(CMRO2n,8,8);


end

time=-1:0.1:4;
imagesc(ATP2d)