%  close all;
clc;
clear;

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
% area_select1=zeros(caps,caps);
% area_select1(31-6:31+6,35-6:35+6)=1; % 800mm2 = 12x12 pixels approx 1 barrel


lp=1;
ATPref=2.5;
loopnr=0;
for j=1:length(nm)
    
    
    i=nm(j);
    loopnr=loopnr+1;
    
    V=x(i,1:nE)';
 
    V1(1,:)=V0;
    V1(2,:)=V;
  
 newd0=2*sqrt((V0*1e9)./(pi*L0));
R(1,:)= newd0/2;

newdn=2*sqrt((V*1e9)./(pi*L0));
R(2,:)=newdn/2;
    HbT=R;

% for trial=1:2
%     if trial==1
%     capsHb2d=caps2D(HbT(end,lev==3),3,0);
%     capsHb2d0=caps2D(HbT(1,lev==3),3,0);
%     else
        capsHbT2d=caps2D(HbT(end,lev==3),3,0);
    capsHbT2d0=caps2D(HbT(1,lev==3),3,0);

%     end
%     area_select1=ones(size(capsHb2d));
    g=g+0.1;
      delHbT= (capsHbT2d-capsHbT2d0)./capsHbT2d0;
      delHbT=round(delHbT,2);
     nonz=sum(delHbT(delHbT>0))/numel(find((delHbT>0)));
%     Avg_Toxy=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
%  
%     Avg_Toxy0=sum(sum( capsHb2d0.*area_select1)/numel(find(area_select1>0)));

    Toxy(loopnr)= 100*nonz;%(Avg_Toxy- Avg_Toxy0)/ Avg_Toxy0;
  
    
   
end

time=-1:0.1:4;


figure(2);plot(time, Toxy);hold on;
% end

title('Change in Radius')
xlabel('TIme in seconds'); ylabel('Percentage change');
