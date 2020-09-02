%  close all;
clc;
clear;
inside=0;
border=0;
outside=0;
loop=1;
for pl=0:0
point=[33+pl,30];
for avg=1:loop
for num=1:1
% filename = 'K1_2cap_2art_K2_0.05_0.5s_cmro2delay_new_pv_reltn_1.7_exponent_0.5s_tau_2Cap_2_art_betadelay.mat';%sprintf('r%d.mat',inps(num)) ;
 filename = 'FINAL_amp3.mat';
 load(filename);k=size(x,1);
      nm=1:k;
load('XY.mat');
load('vessel_charac.mat'); 


neu_length=neu_size(1)*neu_size(2);

g=-0.1;


nC=layers(4);

beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n0=x(1,3*nE+neu_length+nE+1:end);

caps=sqrt(nC);
area_select1=zeros(caps,caps);
%  area_select1(26:36,30:40)=1;
area_select1(31-2:31+2,35-2:35+2)=1; % 300mm2 around centre of PRINCIPAL BARREL
area_select2=zeros(neu_size);
% area_select2(1:20,1:20)=1;
area_select2(17:21,9:13)=1; % FAR AWAY 
area_select3=zeros(neu_size);
% area_select3(25:37,24:29)=1;
area_select3(23:27,35-2:35+2)=1;%BORDER

figure;imagesc(area_select1+area_select2+area_select3)
lp=1;
ATPref=2.5;
loopnr=0;
for j=1:length(nm)
   
  
    i=nm(j);
    loopnr=loopnr+1;

V=x(i,1:nE);
S=x(i,nE+1:2*nE);
PO2tis=x(i,2*nE+1:3*nE);
ATP=x(i,3*nE+1:3*nE+neu_length);
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

 capsHb2d=caps2D(HbT(end,lev==4),4,0);
  capsHb2d0=caps2D(HbT(1,lev==4),4,0);
 
g=g+0.1;

 Avg_Hbt_1=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
Avg_HbT_2=sum(sum(capsHb2d.*area_select2))/numel(find(area_select2>0));
Avg_HbT_3=sum(sum(capsHb2d.*area_select3))/numel(find(area_select3>0));
 Avg_Hbt0_1=sum(sum( capsHb2d0.*area_select1)/numel(find(area_select1>0)));
Avg_HbT0_2=sum(sum(capsHb2d0.*area_select2))/numel(find(area_select2>0));
Avg_HbT0_3=sum(sum(capsHb2d0.*area_select3))/numel(find(area_select3>0));%  %   [delF,delV,delS,delHbT,delHbO,delHb]=percnt_chng(Fin1,V1,S1,HbT,HbO,Hb,X,Y,lev,bndry,j,CMRO2n,BETA,CMRO2n0); pause(0.1);

Avg_Hbtt_1(loopnr)= (Avg_Hbt_1- Avg_Hbt0_1)/ Avg_Hbt0_1;
Avg_HbTt_2(loopnr)=(Avg_HbT_2- Avg_HbT0_2)/ Avg_HbT0_2;
Avg_HbTt_3(loopnr)=(Avg_HbT_3- Avg_HbT0_3)/ Avg_HbT0_3;

end
end
time=-1:0.1:4;
inside=inside+Avg_Hbtt_1;
outside=outside+Avg_HbTt_2;
border=border+Avg_HbTt_3;
end
%  set(groot,'defaultLineLineWidth',2)
figure();plot(time, inside/loop);hold on;
plot( time,border/loop,'g');hold on;
plot( time,outside/loop,'r');
legend('Center of C3 barrel','On the border of C3','Far away from C3')
end

