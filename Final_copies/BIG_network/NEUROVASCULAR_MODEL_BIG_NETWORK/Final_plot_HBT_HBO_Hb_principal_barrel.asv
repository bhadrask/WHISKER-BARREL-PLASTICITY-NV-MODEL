%  close all;
clc;
clear;
load('HbO_exp_v1.csv');
EHbOx=(HbO_exp_v1(:,1));
EHbOy=100*HbO_exp_v1(:,2);
load('HbT_exp.csv');
EHbTx=HbT_exp(:,1);
EHbTy=100*HbT_exp(:,2);
load('Hb_exp.csv');
EHbx=Hb_exp(:,1);
EHby=100*Hb_exp(:,2);

loop=1;

% filename = 'K1_2cap_2art_K2_0.05_0.5s_cmro2delay_new_pv_reltn_1.7_exponent_0.5s_tau_2Cap_2_art_betadelay.mat';%sprintf('r%d.mat',inps(num)) ;
filename = 'final_amp3.mat';
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
area_select1(31-6:31+6,35-6:35+6)=1; % 800mm2 = 12x12 pixels approx 1 barrel


lp=1;
ATPref=2.5;
loopnr=0;
for j=1:length(nm)
    
    
    i=nm(j);
    loopnr=loopnr+1;
    
    V=x(i,1:nE);
    S=x(i,nE+1:2*nE);
    V1=zeros(2,size(Fin0,1));
    S1=zeros(2,size(Fin0,1));
    
    V1(1,:)=V0';
    V1(2,:)=V';
    S1(1,:)=S0';
    S1(2,:)=S';
    HbT=2.3*V1;
    HbO=2.3.*S1.*V1;
    Hb=2.3*(1-S1).*V1;
    
    
    capsHb2d=caps2D(HbT(end,lev==4),4,0);
    capsHb2d0=caps2D(HbT(1,lev==4),4,0);
    oxyH=caps2D(HbO(end,lev==4),4,0);
    oxyH0=caps2D(HbO(1,lev==4),4,0);
    deoxyH=caps2D(Hb(end,lev==4),4,0);
    deoxyH0=caps2D(Hb(1,lev==4),4,0);
    
    g=g+0.1;
    
    Avg_Toxy=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
    Avg_oxy=sum(sum(oxyH.*area_select1))/numel(find(area_select1>0));
    Avg_deoxy=sum(sum(deoxyH.*area_select1))/numel(find(area_select1>0));
    Avg_Toxy0=sum(sum( capsHb2d0.*area_select1)/numel(find(area_select1>0)));
    Avg_oxy0=sum(sum(oxyH0.*area_select1))/numel(find(area_select1>0));
    Avg_deoxy0=sum(sum(deoxyH0.*area_select1))/numel(find(area_select1>0));
    Toxy(loopnr)= 100*(Avg_Toxy- Avg_Toxy0)/ Avg_Toxy0;
    oxy(loopnr)=100*(Avg_oxy- Avg_oxy0)/ Avg_oxy0;
    deoxy(loopnr)=100*(Avg_deoxy- Avg_deoxy0)/ Avg_deoxy0;
    
   
end

time=-1:0.1:4;


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
