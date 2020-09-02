close all;
clc;
clear;
filename = 'FINAL_amp3.mat';
load(filename);k=size(x,1);
input1=IN.input;
nm=1:k;

load('XY.mat');
load('vessel_charac.mat');
neu_length=neu_size(1)*neu_size(2);
eA=bndry(5); eV=bndry(6);
g=-0.1;
nC=layers(4);

beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n0=x(1,3*nE+neu_length+nE+1:end);

caps=sqrt(nC);
area_select1=zeros(caps,caps);
area_select1(31-6:31+6,35-6:35+6)=1;
neu_select1=zeros(neu_size);
neu_select1(29:33,33:37)=1;
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
    
    beta=x(i,3*nE+neu_length+1:3*nE+neu_length+nE);
    CMRO2n=x(i,3*nE+neu_length+nE+1:end);
    dVdt=zeros(size(V));
    PO2tissn1=x(i,2*nE+1:3*nE);
    PO2tissn=PO2tissn1(eA+1:eV-1);
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
    
    if (g<1) || (g>2)
        input=zeros(size(input1));
    else
        input=input1;
    end
    
    Fpoint=caps2D(F(end,lev==4),4,0);
    Vpoint=caps2D(V1(end,lev==4),4,0);
    Fpoint0=caps2D(F(1,lev==4),4,0);
    Vpoint0=caps2D(V1(1,lev==4),4,0);
    g=g+0.1;
    ATP2d=reshape(ATP,64,64);
    CMRO2_2d=reshape(CMRO2n,64,64);
    avg_CMRO2n=sum(sum(CMRO2_2d.*area_select1)/numel(find(area_select1>0)));
    CMRO20_2d=reshape(CMRO2n0,64,64);
    avg_CMRO2n0=sum(sum(CMRO20_2d.*neu_select1)/numel(find(neu_select1>0)));
    AVG_ATP=sum(sum( ATP2d.*neu_select1)/numel(find(neu_select1>0)));
    ATPper(loopnr)=100*((AVG_ATP)-2.5)/2.5;
    CMR02_per(loopnr)=avg_CMRO2n/avg_CMRO2n0;%100*(avg_CMRO2n-avg_CMRO2n0)/avg_CMRO2n0;
    Fp=sum(sum(Fpoint.*area_select1)/numel(find(area_select1>0)));
    Fp0=sum(sum(Fpoint0.*area_select1)/numel(find(area_select1>0)));
    Flow(loopnr)=100*(Fp- Fp0)/ Fp0;
    Vp=sum(sum(Vpoint.*area_select1)/numel(find(area_select1>0)));
    Vp0=sum(sum(Vpoint0.*area_select1)/numel(find(area_select1>0)));
    Vol(loopnr)=100*(Vp- Vp0)/ Vp0;
    
end

time=-1:0.1:4;


figure();plot(time, Flow,'r');
hold on;
plot( time,Vol,'b');
hold on;

legend('CBF','CBV');
title('CBF and CBV')
xlabel('Time in seconds'); ylabel('Percentage change')


for i=1:size(neuout,3)
    neunet(i)=sum(sum( neuout(:,:,i).*area_select1)/numel(find(area_select1>0)));
end
figure();stem(Ttime,neunet);title('Neural Output');
figure();plot(time,ATPper);

xlabel('Time in seconds');ylabel('% drop in ATP')
title('Percentage drop in concentration of ATP')
figure();plot(time,CMR02_per);title('Percentage change in concentration of CMRO2')