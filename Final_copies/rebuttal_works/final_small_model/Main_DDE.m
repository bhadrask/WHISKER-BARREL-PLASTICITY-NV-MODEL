clc
clear all
close all
tic;
%% Load all the saved mat files going to be used
load('vessel_charac.mat');                  % This loads the diameter, length and the distribution of arteries, veins and capilliaries
load('XY.mat');
load('lissom_10atpmax_1p_6q_6r.mat');                            % this describes the location of vessel branches in a 2D matrix
file_vessel='XY.mat';
file_vessel_chara='vessel_charac.mat';
load X_test_4vasc_amp5.mat;                      % This loads the input
opts.niter=1; % redefine niter since the network is already trained
opts.etaaff=[0.01];
opts.etaexc=[0.005];
opts.etainhib=[0.005];
%% Initialize the vessel parameters
nE=length(E); nN=N2(end)+1; eA=bndry(5); eV=bndry(6);
caps=eV-eA-1;
 d0=0.85*d0;
d2=d0(1);
% d0(1)=0.4*d0(1);
R0=128*15*(L0/2)*1e3./(pi*d0.^4);%units in mmHg/microlitres
V0=1e-9*pi*L0.*d0.^2/4; %units in microlitres
%% Set up the neural layer with pretrained weights

Amax=ATPmax;
% setGlobalx(lissom.layers{2}.Zold,lissom.layers{2}.Zold);

neu_size=lissom.layers{2}.dim;
% ATPref=2.5;
ATP=ATPref* ones(lissom.layers{2}.dim(1)* lissom.layers{2}.dim(2),1);
thr_atp=(ATPref-(0.10*ATPref));
%% Choose an input

input=X1(:,:,1:3);%(:,:,1);
%% Initializing vessel parameters
[Ae,Pn0,Pe0,Fin0,Fout0]=vsc_func0a(E,N1,N2,bndry,R0,V0,para); % initializing vascular parameter, pressure and flow
PO2tis=15*ones(nE,1);
PO2tismax=PO2tis;
PO2tis4art=30;
W_gauss=gaussian_filter_wt(3,1);
loc_vas=findloc(file_vessel,file_vessel_chara);
PO2tisv=project2vasc(bndry,PO2tis,PO2tis4art,W_gauss,loc_vas);

Sen=0.94;
%  load('S0_gauss_4,1.mat');


options1 = optimoptions('fsolve','Algorithm','levenberg-marquardt');
% [S0,fval]=fsolve(@(S) vsc_func0b_Avg_aditi_new(S,E,N1,N2,Fin0,Fout0,L0,d0,PO2tisv,Sen),linspace(0.9,0.4,nE)',options1);
% S0(S0<0)=0;
load('S0_nov_8_2020','S0');
load('OE_healthy');
OE0=OEhealthy;

[OE,OEn,CMRO2n]=init_func_mod_final(S0,bndry,neu_size,PO2tisv,L0,d0,W_gauss,file_vessel,file_vessel_chara,OE0);
% OEhealthy=OEn;
% save('OE_healthy','OEhealthy');
% OE0=OEn;
beta0=Ae;
CMRO2nref=CMRO2n;
x0=[V0;S0;PO2tis;ATP;beta0;CMRO2n];


%% the ODE function
% options2 = odeset('RelTol',1e-6,'AbsTol',1e-6);
% [t,x]=ode45(@(t,x)  function4ode(t,x,Ae,lissom,opts,input,CMRO2nref,ATPref,Z_constraint,PO2tis4art,W_gauss,loc_vas,Sen,file_vessel,OE0,file_vessel_chara),0:0.1:5,x0,options2);
training=1;


% cc1=0;cc2=0;
% setGlobalx(cc1,cc2);
neu_old=zeros(neu_size);
c=0;
%  setGlobalneu(neu_old,c);
tau=0.5;
lags=tau*ones(size(x0));
% lags(3*nE+1:3*nE+neu_length)=1;
ckc=1;
neu_size=lissom.layers{2}.dim;
neu_length=neu_size(1)*neu_size(2);
% Waffref=cell2mat(lissom.layers{2}.waff);
for epoch=1:4
    epoch
    pos=randperm(size(input,3));
    
    
    for jj=1:length(pos)
        lissom.layers{2}.Zold =zeros(lissom.layers{2}.dim(1)*lissom.layers{2}.dim(2),1);
        lissom.layers{2}.Z  =zeros(lissom.layers{2}.dim(1)*lissom.layers{2}.dim(2),1);
        lissom.layers{2}.Zprev  =zeros(lissom.layers{2}.dim(1)*lissom.layers{2}.dim(2),1);
        setGlobalx(lissom,training);
        setGlobaly(lissom);
        neuout=zeros(lissom.layers{2}.dim(1)*lissom.layers{2}.dim(2),length(0:0.1:4));
        cnt=1;
        neutime=zeros(1,size(neuout,2));
        setGlobalneu(cnt,neuout,neutime);
        input1=input(:,:,pos(jj));
        tend=2;
        steps=0.1;
        time_span=0:steps:tend;
        options2 = ddeset('RelTol',1e-4,'AbsTol',1e-4);
        options3 = odeset(options2,'NonNegative',1);
        [x]=dde23(@(t,x,z)  function4DDE(t,x,z,Ae,opts,input1,CMRO2nref,ATPref,PO2tis4art,W_gauss,loc_vas,Sen,file_vessel,file_vessel_chara,training,lissom,tau,Fin0,Fout0,V0,PO2tismax,ATPmax,thr_atp),lags,x0,time_span,options3);
        tint = linspace(0,tend,length(time_span));
        Sint = deval(x,tint);
        
        x0=Sint(:,end);
        
        lissom=getGlobaly;
        
        [cnt,neuout,Ttime]=getGlobalneu;
        x2=Sint';
        ATP=x2(end,3*nE+1:3*nE+neu_length);
        MinATP(ckc)=min(ATP);
        ckc=ckc+1;
    end
end
x=Sint';
save('TEST1.mat','lissom','ATPmax','x','E','N1','N2','bndry','para','neu_size','nE','L0','S0','Fin0','Fout0','V0','cnt','neuout','Ttime','MinATP','thr_atp','opts');
% plot4mx_imagesc(x,E,N1,N2,bndry,para,neu_size,nE,L0,S0,Fin0,Fout0,V0,input1)
yo=toc;
if yo>60 && yo<=3600
    yt=yo/60;
    ds=['elapsed time = ',num2str(yt),'min'];
    disp(ds);
elseif yo>3600
    yt=yo/(60*60);
    ds=['elapsed time = ',num2str(yt),'hours'];
    disp(ds);
else
    ds=['elapsed time = ',num2str(yo),'seconds'];
    disp(ds);
end
% load('X_test_4vasc.mat');
Xtest=X1;
neu_size=lissom.layers{2}.dim;
neu_length=neu_size(1)*neu_size(2);
ATP=x(end,3*nE+1:3*nE+neu_length);
min(ATP)
% [yy,mean_act]=inst_response(lissom,opts,Xtest,label,ATP,ATPmax,thr_atp);
MAPS(lissom,opts,Xtest,label,ATP,ATPmax,thr_atp);
mean_act