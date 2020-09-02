clc;
clear;
close all;
tic;
%% Load all the saved mat files going to be used
load('vessel_charac.mat');                  % This loads the diameter, length and the distribution of arteries, veins and capilliaries
load('XY.mat');                             % this describes the location of vessel branches in a 2D matrix
file_vessel='XY.mat';
file_vessel_chara='vessel_charac.mat';
load Z_constraint.mat;                      % This describes the shape of the neural layer
load('amp14_X_test_sigma_4.mat');                      % This loads the input
load('lissom_10atpmax_newopts_14q_8r.mat');                % Loads the trained weights
opts.niter=1; % redefine niter since the network is already trained
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
ATPref=2.5;
ATP=ATPref*ones(lissom.layers{2}.dim(1)* lissom.layers{2}.dim(2),1);

%% Input stimulus definition
i1=12;                                       % choose the input. put in a loop if needed
IN.input=X1(:,:,i1);
IN.Tstart=1;  % time at which input stimulus is given
IN.Tstop=2; % time at which input is taken off
IN.Tau=0.5; % time taken for Neural signal to cause vasodilation

%% Initializing vessel parameters
[Ae,Pn0,Pe0,Fin0,Fout0]=vsc_func0a(E,N1,N2,bndry,R0,V0,para); % initializing vascular parameter, pressure and flow
PO2tis=15*ones(nE,1);
PO2tismax=PO2tis;
PO2tis4art=30;
W_gauss=gaussian_filter_wt(10,2);
loc_vas=findloc(file_vessel);
PO2tisv=project2vasc(bndry,PO2tis,PO2tis4art,W_gauss,loc_vas);
Sen=0.94;
%% steps to find initial saturation
% since this is time consuming, we train it once and then load the saved mat file
% options1 = optimoptions('fsolve','Algorithm','levenberg-marquardt');
% [S0,fval]=fsolve(@(S) vsc_func0b_Avg_aditi_new(S,E,N1,N2,Fin0,Fout0,L0,d0,PO2tisv,Sen),linspace(0.9,0.4,nE)',options1);
% save('S0_gauss_10_2','S0');
load('S0_gauss_10_2.mat');
%% Other vessel parameters
load('OE_healthy');
OE0=OEhealthy;
[OE,OEn,CMRO2n]=init_func_mod_final_mod(S0,bndry,neu_size,PO2tisv,L0,d0,W_gauss,file_vessel,OE0);
beta0=Ae;
CMRO2nref=CMRO2n;
x0=[V0;S0;PO2tis;ATP;beta0;CMRO2n];

%% the ODE function
options2 = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,x]=ode45(@(t,x)  function4ode(t,x,Ae,lissom,opts,Amax,IN,CMRO2nref,ATPref,Z_constraint,PO2tis4art,W_gauss,loc_vas,Sen,file_vessel,file_vessel_chara,PO2tismax),0:0.1:5,x0,options2);
[cnt,neuout,Ttime]=getGlobalneu;
save('FINAL_amp14.mat','x','E','N1','N2','bndry','para','neu_size','nE','L0','S0','Fin0','Fout0','V0','IN','cnt','neuout','Ttime');
%worked_lissom_10atpmax_newopts_14q_8r__c1_10_t1_0.5_t2_0.1_tau_0.5_c3_0.7_t3_2_tsat_1.mat
%% Estimating computational time
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