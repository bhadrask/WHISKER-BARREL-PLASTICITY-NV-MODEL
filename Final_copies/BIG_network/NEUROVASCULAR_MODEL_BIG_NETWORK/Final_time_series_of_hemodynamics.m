% close all;
clc;
clear all;

filename = 'final_amp25.mat';
load(filename);k=size(x,1);
load('XY.mat');
load('vessel_charac.mat');

input1=IN.input;

nm=1:k;

neu_length=neu_size(1)*neu_size(2);
eA=bndry(5); eV=bndry(6);

g=-200;
beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n0=x(1,3*nE+neu_length+nE+1:end);
lp=1;
ATPref=2.5;
loopnr=0;
 vec=9:2:length(nm);
    for j=vec
    
    
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
    if (g<0) || (g>1000)
        input=zeros(size(input1));
    else
        input=input1;
    end
    capsHbT2d=caps2D(HbT(end,lev==4),4,0);
    capsHbT2d0=caps2D(HbT(1,lev==4),4,0);
    capsHbO2d=caps2D(HbO(end,lev==4),4,0);
    capsHbO2d0=caps2D(HbO(1,lev==4),4,0);
    capsHb2d=caps2D(Hb(end,lev==4),4,0);
    capsHb2d0=caps2D(Hb(1,lev==4),4,0);
    %% Calculating and Plotting fractional change
    delHbT= (capsHbT2d-capsHbT2d0)./capsHbT2d0;
    delHbO= (capsHbO2d-capsHbO2d0)./capsHbO2d0;
    delHb= (capsHb2d-capsHb2d0)./capsHb2d0;
    index=find(lev==4);
         figure(4);
        ax(loopnr)= subplot(6,5,loopnr);
        delHbT=round(delHbT,2);
        delHbO=round(delHbO,2);
        delHb=round(delHb,2);
        imagesc(delHbT);axis off; 
        h=colorbar;caxis([-0.15 0.151]);
        set(h, 'Position', [.9 .11 .02 .8150]); %[left, bottom, width, height]
        pos=get(ax(loopnr), 'Position');
        set(ax(loopnr), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
        title([num2str(g/1000),'s']);
       
        figure(5);  ax(loopnr)=subplot(6,5,loopnr);
        imagesc(delHbO); axis off;
         h=colorbar;caxis([-0.15 0.151]);
        set(h, 'Position', [.9 .11 .02 .8150]); %[left, bottom, width, height]
        pos=get(ax(loopnr), 'Position');
        set(ax(loopnr), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
        title([num2str(g/1000),'s']);
        
        figure(6);  ax(loopnr)=subplot(6,5,loopnr);
        imagesc(delHb);axis off;
         h=colorbar;caxis([-0.15 0.151]);
        set(h, 'Position', [.9 .11 .02 .8150]); %[left, bottom, width, height]
        pos=get(ax(loopnr), 'Position');
        set(ax(loopnr), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
        title([num2str(g/1000),'s']);
        if g==1600
           figure(1); imagesc(delHbO)
        end
        g=g+200;
    %   [delF,delV,delS,delHbT,delHbO,delHb]=percnt_chng(Fin1,V1,S1,HbT,HbO,Hb,X,Y,lev,bndry,j,CMRO2n,BETA,CMRO2n0); pause(0.1);
    
end

figure(4);suptitle('HbT');
figure(5);suptitle('HbO');
figure(6);suptitle('Hb');

