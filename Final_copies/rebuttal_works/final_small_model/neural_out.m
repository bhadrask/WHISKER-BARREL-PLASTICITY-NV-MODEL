clc;
clear;
close all;
tic;
%% Load all the saved mat files going to be used
load('vessel_charac.mat');                  % This loads the diameter, length and the distribution of arteries, veins and capilliaries
load('XY.mat'); 
load('lissom_10atpmax_1p_10q_8r');                            % this describes the location of vessel branches in a 2D matrix
file_vessel='XY.mat';  
file_vessel_chara='vessel_charac.mat';
load('X_test_4vasc.mat');                      % This loads the input


opts.niter=1; % redefine niter since the network is already trained
 % redefine niter since the network is already trained
% opts.p=2; opts.q=10; opts.r=5;
%% Initialize the vessel parameters
nE=length(E); nN=N2(end)+1; eA=bndry(5); eV=bndry(6);
caps=eV-eA-1;
d0=0.85*d0;
R0=128*15*(L0/2)*1e3./(pi*d0.^4);%units in mmHg/microlitres
V0=1e-9*pi*L0.*d0.^2/4; %units in microlitres
%% Set up the neural layer with pretrained weights
lissom.layers{2}.Zold =zeros(lissom.layers{2}.dim);
lissom.layers{2}.Z  =zeros(lissom.layers{2}.dim);
Amax=ATPmax;
setGlobalx(lissom.layers{2}.Zold,lissom.layers{2}.Zold);
neuout=zeros(lissom.layers{2}.dim(1),lissom.layers{2}.dim(2),length(0:0.1:4));
cnt=1;
neutime=zeros(1,size(neuout,3));

setGlobalneu(cnt,neuout,neutime);
neu_size=lissom.layers{2}.dim;
Z_constraint=ones(neu_size);
ATPref=2.5;
ATP=ATPref*ones(lissom.layers{2}.dim(1)* lissom.layers{2}.dim(2),1);

%% Input stimulus definition
                                       % choose the input. put in a loop if needed
IN.input=X1(:,:,1);
IN.Tstart=1;  % time at which input stimulus is given
IN.Tstop=2; % time at which input is taken off
IN.Tau=0.8; % time taken for Neural signal to cause vasodilation

%% Initializing vessel parameters
[Ae,Pn0,Pe0,Fin0,Fout0]=vsc_func0a(E,N1,N2,bndry,R0,V0,para); % initializing vascular parameter, pressure and flow
PO2tis=15*ones(nE,1);
PO2tismax=PO2tis;
PO2tis4art=30;
W_gauss=gaussian_filter_wt(3,1);
loc_vas=findloc(file_vessel);
%% Other vessel parameters

beta0=Ae;

x0=beta0;

%% the ODE function
options2 = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,x]=ode45(@(t,x)  fun_neu_out(t,x,Ae,lissom,opts,Amax,IN,Z_constraint,W_gauss,loc_vas,file_vessel,file_vessel_chara),0:0.1:5,x0,options2);  
 %% time
 yo=toc;
if yo>60 && yo<=3600
    yt=yo/60;
    ds=['elapsed time = ',num2str(yt),'min'];
    disp(ds);
elseif yo>3600
    yt=yo/(60*60);
    ds=['elapsed time = ',num2str(yt),'hours'];
    disp(ds);
end

%% plot
 nC=layers(4);
 caps=sqrt(nC);
area_select1=zeros(caps,caps);
area_select1(1:4,1:4)=1;
 beta0=caps2D(x(1,lev==4),4,0);
for i=1:size(x,1)
beta=caps2D(x(i,lev==4),4,0);

 Bt=sum(sum( beta.*area_select1)/numel(find(area_select1>0)));
  Bt0=sum(sum( beta0.*area_select1)/numel(find(area_select1>0)));
  perbta(i)=100*(Bt-Bt0)/Bt0;
 
end
[cnt,neuout,Ttime]=getGlobalneu;
for i=1:size(neuout,3)
%     NET(:,:,i)=conv2(neuout(:,:,i),W_gauss,'same');
     neunet(i)=sum(sum( neuout(:,:,i).*area_select1)/numel(find(area_select1>0)));
end
figure();plot(0:0.1:5,perbta);title('beta');
figure();stem(neunet);title('neu');