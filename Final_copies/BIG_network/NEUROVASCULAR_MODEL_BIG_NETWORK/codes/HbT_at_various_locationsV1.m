close all;
clc;
clear all;
% load('output_4_plot_c3_amp=100.mat');
% disp('load the saved mat file');
% keyboard; 
inps=800;%[10,50,80,100,200,300,400,500,800];
inside=0;
border=0;
outside=0;
loop=1;
for pl=0:0
point=[33+pl,30];
for avg=1:loop
for num=1:1
filename = 'FINAL_amp15.mat';% sprintf('r%d.mat',inps(num)) ;
load(filename);k=size(x,1);
load('XY.mat');
load('vessel_charac.mat'); 
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


nC=layers(4);

caps=sqrt(nC);
area_select1=zeros(caps,caps);
area_select1(26:36,30:40)=1;
area_select2=zeros(neu_size);
area_select2(1:20,1:20)=1;
area_select3=zeros(neu_size);
area_select3(25:37,24:29)=1;

beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n0=x(1,3*nE+neu_length+nE+1:end);
lp=1;
ATPref=2.5;
loopnr=0;
for j=1:length(nm)
   
  
    i=nm(j);
    loopnr=loopnr+1;
V=x(i,1:nE)';
newd=2*sqrt((V*1e9)./(pi*L0));
 R=128*15*(L0/2)*1e3./(pi*newd.^4);

S=x(i,nE+1:2*nE)';
ATP=x(i,3*nE+1:3*nE+neu_length);
% beta=x(i,3*nE+neu_length+1:end);
beta=x(i,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n=x(i,3*nE+neu_length+nE+1:end);

V1=zeros(2,size(Fin0,1));
S1=zeros(2,size(Fin0,1));

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

%  capsHb2d=reshape(CMR(2,:),8,8);
%   capsHb2d0=reshape(CMR(1,:),8,8);
%% Calculating and Plotting fractional change
% delHbT= 100*(HbT-repmat(HbT(1,:),size(V1,1),1))./repmat(HbT(1,:),size(V1,1),1);
% delHbO= 100*(HbO-repmat(HbO(1,:),size(V1,1),1))./repmat(HbO(1,:),size(V1,1),1);
% delHb= 100*(Hb-repmat(Hb(1,:),size(V1,1),1))./repmat(Hb(1,:),size(V1,1),1);
% delS= 100*(S1-repmat(S1(1,:),size(S1,1),1))./repmat(S1(1,:),size(S1,1),1);

% caps_final=capsHb2d.*area_select1;%./max(max(capsHb2d.*area_select));
%  Avg_Hbt_1(loopnr)=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
% Avg_HbT_2(loopnr)=sum(sum(capsHb2d.*area_select2))/numel(find(area_select2>0));
% Avg_HbT_3(loopnr)=sum(sum(capsHb2d.*area_select3))/numel(find(area_select3>0));
 Avg_Hbt_1=sum(sum( capsHb2d.*area_select1)/numel(find(area_select1>0)));
Avg_HbT_2=sum(sum(capsHb2d.*area_select2))/numel(find(area_select2>0));
Avg_HbT_3=sum(sum(capsHb2d.*area_select3))/numel(find(area_select3>0));
 Avg_Hbt0_1=sum(sum( capsHb2d0.*area_select1)/numel(find(area_select1>0)));
Avg_HbT0_2=sum(sum(capsHb2d0.*area_select2))/numel(find(area_select2>0));
Avg_HbT0_3=sum(sum(capsHb2d0.*area_select3))/numel(find(area_select3>0));%  %   [delF,delV,delS,delHbT,delHbO,delHb]=percnt_chng(Fin1,V1,S1,HbT,HbO,Hb,X,Y,lev,bndry,j,CMRO2n,BETA,CMRO2n0); pause(0.1);
 meanatp(loopnr)=100*(mean(ATP)-mean(ATPref))/mean(ATPref);
  meanatp1(loopnr)=mean(ATP);
Avg_Hbtt_1(loopnr)= 100*(Avg_Hbt_1- Avg_Hbt0_1)/ Avg_Hbt0_1;
Avg_HbTt_2(loopnr)=100*(Avg_HbT_2- Avg_HbT0_2)/ Avg_HbT0_2;
Avg_HbTt_3(loopnr)=100*(Avg_HbT_3- Avg_HbT0_3)/ Avg_HbT0_3;
 Toxy(loopnr)=100*(capsHb2d(point(1),point(2))-capsHb2d0(point(1),point(2)))/capsHb2d0(point(1),point(2));
 oxy(loopnr)=100*(oxyH(point(1),point(2))-oxyH0(point(1),point(2)))/oxyH0(point(1),point(2));
 deoxy(loopnr)=100*(deoxyH(point(1),point(2))-deoxyH0(point(1),point(2)))/deoxyH0(point(1),point(2));
 one(loopnr)=capsHb2d(point(1),point(2));
 two(loopnr)=oxyH(point(1),point(2));
 three(loopnr)=deoxyH(point(1),point(2));
 


end
end
time=-1:0.1:4;
% plot(time, Avg_Hbtt_1);hold on;
% plot( time,Avg_HbTt_2,'r');hold on;
% plot( time,Avg_HbTt_3,'g');
inside=inside+Avg_Hbtt_1;
outside=outside+Avg_HbTt_2;
border=border+Avg_HbTt_3;
end
%  figure();plot(time,meanatp);
%   figure();plot(time,meanatp1);
figure(1);plot(time, inside/loop);hold on;
plot( time,border/loop,'g');hold on;
plot( time,outside/loop,'r');
legend('Center of C3 barrel','On the border of C3','Far away from C3')
figure();plot(time, Toxy,'r');hold on;
plot( time,oxy,'b');hold on;
plot( time,deoxy,'g');
legend('HbT','HbO','Hb');
title('Hemoglobin oxygenation profile at the centre of the barrel')
xlabel('TIme in seconds'); ylabel('Percentage change')

end


