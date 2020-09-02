  close all;
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
 filename = 'FINAL_amp8.mat';
 load(filename);k=size(x,1);
load('XY.mat');
load('vessel_charac.mat');  
input1=IN.input;
     nm=1:k;

neu_length=neu_size(1)*neu_size(2);
eA=bndry(5); eV=bndry(6);

g=-0.1;

nC=layers(4);

beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n0=x(1,3*nE+neu_length+nE+1:end);

caps=sqrt(nC);
area_select1=zeros(caps,caps);
  area_select1(24:40,25:32)=1;
%  area_select1(25:37,29:41)=1;

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

g=g+0.1;


 Toxy_per=100*((capsHb2d-capsHb2d0)./capsHb2d0);
 oxy_per=100*((oxyH-oxyH0)./oxyH0);
 deoxy_per=100*((deoxyH-deoxyH0)./deoxyH0);
 Toxy_per(abs(Toxy_per)<1)=0;
 oxy_per(abs(oxy_per)<1)=0;
  deoxy_per(abs(deoxy_per)<1)=0;
if round(g,1)==2.6
figure;subplot 311;imagesc(Toxy_per);
subplot 312;imagesc(oxy_per);
subplot 313;imagesc(deoxy_per);
end
 Toxy(loopnr)=sum(sum(Toxy_per))/(30*30);% sum(sum(Toxy_per.*area_select1)/numel(find(area_select1>0)));
 oxy(loopnr)= sum(sum(oxy_per))/(30*30);%sum(sum(oxy_per.*area_select1)/numel(find(area_select1>0)));
 deoxy(loopnr)=sum(sum(deoxy_per))/(30*30);%sum(sum(deoxy_per.*area_select1)/numel(find(area_select1>0)));
 

end
end
time=-1:0.1:4;

end
%  set(groot,'defaultLineLineWidth',3.0)
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

% figure();plot(time, Flow,'r');
% hold on;
% plot( time,Vol,'b');
% hold on;
% 
% legend('CBF','CBV');
% title('CBF and CBV')
% xlabel('Time in seconds'); ylabel('Percentage change')
end

for i=1:size(neuout,3)
%     NET(:,:,i)=conv2(neuout(:,:,i),W_gauss,'same');
     neunet(i)=sum(sum( neuout(:,:,i).*area_select1)/numel(find(area_select1>0)));
end
% figure();stem(Ttime,neunet);title('neu');
% figure();plot(time,ATPper); 
% % ylim([2,2.9]) ; 
% xlabel('Time in seconds');ylabel('% drop in ATP')
% title('Percentage drop in concentration of ATP')
% figure();plot(time,CMR02_per);