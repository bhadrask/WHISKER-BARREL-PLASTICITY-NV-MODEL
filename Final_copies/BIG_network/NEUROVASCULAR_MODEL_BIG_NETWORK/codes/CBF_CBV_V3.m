%  close all;
clc;
clear all;
load('HbO_exp_v1.csv');
EHbOx=(HbO_exp_v1(:,1));
EHbOy=100*HbO_exp_v1(:,2);
load('HbT_exp.csv');
EHbTx=HbT_exp(:,1);
EHbTy=100*HbT_exp(:,2);
load('Hb_exp.csv');
EHbx=Hb_exp(:,1);
EHby=100*Hb_exp(:,2);
inside=0;
border=0;
outside=0;
loop=1;
for pl=0:0
point=[33+pl,30];
for avg=1:loop
for num=1:1
% filename = 'K1_2cap_2art_K2_0.05_0.5s_cmro2delay_new_pv_reltn_1.7_exponent_0.5s_tau_2Cap_2_art_betadelay.mat';%sprintf('r%d.mat',inps(num)) ;
 filename = 'final_amp25.mat';
 load(filename);k=size(x,1);
load('XY.mat');
load('vessel_charac.mat'); 
% load('lissomtrain1.mat');
% % load('Z_constraint.mat');
% file_neural='lissomtrain1.mat'; 
% % [neu_size]=(size(Z_constraint));
% [lissom,opts] = define(neu_size);
% 
% [lissom] = initialise4usage(lissom,file_neural); 
load('lissom_10atpmax_newopts_14q_8r.mat');  
input1=IN.input;
% if k>20
% nm=1:5:k;
% else
     nm=1:k;
% end
load('XY.mat');
load('vessel_charac.mat'); 
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
area_select1(30:40,26:36)=1;
area_select2=zeros(neu_size);
area_select2(1:20,1:20)=1;
area_select3=zeros(neu_size);
area_select3(1:20,20:50)=1;


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
if (g<1) || (g>2)
    input=zeros(size(input1));
else
    input=input1;
end
 capsHb2d=caps2D(HbT(end,find(lev==4)),4,0);
  capsHb2d0=caps2D(HbT(1,find(lev==4)),4,0);
   oxyH=caps2D(HbO(end,find(lev==4)),4,0);
 oxyH0=caps2D(HbO(1,find(lev==4)),4,0);
  deoxyH=caps2D(Hb(end,find(lev==4)),4,0);
 deoxyH0=caps2D(Hb(1,find(lev==4)),4,0);
 Fpoint=caps2D(F(end,find(lev==4)),4,0);
 Vpoint=caps2D(V1(end,find(lev==4)),4,0);
  Fpoint0=caps2D(F(1,find(lev==4)),4,0);
 Vpoint0=caps2D(V1(1,find(lev==4)),4,0);
g=g+0.1;


 Avg_Hbt_1=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
Avg_HbT_2=sum(sum(capsHb2d.*area_select2))/numel(find(area_select2>0));
Avg_HbT_3=sum(sum(capsHb2d.*area_select3))/numel(find(area_select3>0));
 Avg_Hbt0_1=sum(sum( capsHb2d0.*area_select1)/numel(find(area_select1>0)));
Avg_HbT0_2=sum(sum(capsHb2d0.*area_select2))/numel(find(area_select2>0));
Avg_HbT0_3=sum(sum(capsHb2d0.*area_select3))/numel(find(area_select3>0));%  %   [delF,delV,delS,delHbT,delHbO,delHb]=percnt_chng(Fin1,V1,S1,HbT,HbO,Hb,X,Y,lev,bndry,j,CMRO2n,BETA,CMRO2n0); pause(0.1);

Avg_Hbtt_1(loopnr)= 100*(Avg_Hbt_1- Avg_Hbt0_1)/ Avg_Hbt0_1;
Avg_HbTt_2(loopnr)=100*(Avg_HbT_2- Avg_HbT0_2)/ Avg_HbT0_2;
Avg_HbTt_3(loopnr)=100*(Avg_HbT_3- Avg_HbT0_3)/ Avg_HbT0_3;

 Avg_Toxy=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
Avg_oxy=sum(sum(oxyH.*area_select1))/numel(find(area_select1>0));
Avg_deoxy=sum(sum(deoxyH.*area_select1))/numel(find(area_select1>0));
 Avg_Toxy0=sum(sum( capsHb2d0.*area_select1)/numel(find(area_select1>0)));
Avg_oxy0=sum(sum(oxyH0.*area_select1))/numel(find(area_select1>0));
Avg_deoxy0=sum(sum(deoxyH0.*area_select1))/numel(find(area_select1>0));
 Toxy(loopnr)= 100*(Avg_Toxy- Avg_Toxy0)/ Avg_Toxy0;
 oxy(loopnr)=100*(Avg_oxy- Avg_oxy0)/ Avg_oxy0;
 deoxy(loopnr)=100*(Avg_deoxy- Avg_deoxy0)/ Avg_deoxy0;
 
 Fp=sum(sum(Fpoint.*area_select1)/numel(find(area_select1>0)));
 Fp0=sum(sum(Fpoint0.*area_select1)/numel(find(area_select1>0)));
  Flow(loopnr)=100*(Fp- Fp0)/ Fp0;
   Vp=sum(sum(Vpoint.*area_select1)/numel(find(area_select1>0)));
 Vp0=sum(sum(Vpoint0.*area_select1)/numel(find(area_select1>0)));
 Vol(loopnr)=100*(Vp- Vp0)/ Vp0;

end
end
time=-1:0.1:4;
inside=inside+Avg_Hbtt_1;
outside=outside+Avg_HbTt_2;
border=border+Avg_HbTt_3;
end

figure();plot(time, Toxy,'r');hold on;
plot( time,oxy,'b');hold on;
plot( time,deoxy,'g');
legend('HbT','HbO','Hb');
hold on;
plot(EHbOx,EHbOy,'b.');hold on
plot(EHbTx,EHbTy,'r.');hold on
plot(EHbx,EHby,'g.');
title('Hemoglobin oxygenation profile at the centre of the barrel')
xlabel('TIme in seconds'); ylabel('Percentage change');
% hold on;

figure();plot(time, Flow,'r');
hold on;
plot( time,Vol,'b');
hold on;

legend('CBF','CBV');
title('CBF and CBV')
xlabel('Time in seconds'); ylabel('Percentage change')
end

for i=1:size(neuout,3)
%     NET(:,:,i)=conv2(neuout(:,:,i),W_gauss,'same');
     neunet(i)=sum(sum( neuout(:,:,i).*area_select1)/numel(find(area_select1>0)));
end
figure();stem(Ttime,neunet);title('neu');
