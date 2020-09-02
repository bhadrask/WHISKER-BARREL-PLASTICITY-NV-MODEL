close all;
clc;
clear all;
% load('output_4_plot_c3_amp=100.mat');
% disp('load the saved mat file');
% keyboard;
inps=25;%[10,50,80,100,200,300,400,500,800];

for num=1:1
    % filename = sprintf('r%d.mat',inps(num)) ;
    filename = 'Nat0.5_amp5_tatp2.mat';
    load(filename);k=size(x,1);
    load('XY.mat');
    load('vessel_charac.mat');
    load('lissomtrain1.mat');
    load('Z_constraint.mat');
    file_neural='lissomtrain1.mat';
    [neu_size]=(size(Z_constraint));
    [lissom,opts] = define(neu_size);
    
    [lissom] = initialise4usage(lissom,file_neural);
    input1=input;
    % if k>20
    % nm=1:5:k;
    % else
    nm=1:k;
    % end
    load('XY.mat');
    load('vessel_charac.mat');
    neu_length=neu_size(1)*neu_size(2);
    eA=bndry(5); eV=bndry(6);
    
    g=-200;
    
    nA=sum(layers(1:3));
    nC=layers(4);
    nV=sum(layers(5:7));
    beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nA+nC);
    CMRO2n0=x(1,3*nE+neu_length+nA+nC+1:end);
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
        
        S=x(i,nE+1:2*nE)';
        ATP=x(i,3*nE+1:3*nE+neu_length);
        % beta=x(i,3*nE+neu_length+1:end);
        beta=x(i,3*nE+neu_length+1:3*nE+neu_length+nA+nC);
        CMRO2n=x(i,3*nE+neu_length+nA+nC+1:end);
        dVdt=zeros(size(V));
        PO2tissn1=x(i,2*nE+1:3*nE);
        PO2tissn=PO2tissn1(eA+1:eV-1);
        % [Pn, Pe, Fin, Fout]= vsc_func6(t,R,V,E,N1,N2,bndry,para,beta); %Finding Flow rate
        
        [Pnv,Pe]=vsc_func2_old(E,N1,N2,bndry,R,V,para,beta);
        [PnA,Pe]=vsc_func3(E,N1,N2,bndry,R,Pe,dVdt,para);
        Pn=[PnA;Pnv];
        [dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe);
        Fin1=zeros(2,size(Fin,1));
        Fout1=zeros(2,size(Fin,1));
        Fin1(1,:)=Fin0'; Fout1(1,:)=Fout0';
        Fin1(2,:)=Fin';Fout1(2,:)=Fout';
        F=(Fin1+Fout1)/2;
        V1=zeros(2,size(Fin0,1));
        S1=zeros(2,size(Fin0,1));
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
        capsHbT2d=caps2D(HbT(end,find(lev==4)),4,0);
        capsHbT2d0=caps2D(HbT(1,find(lev==4)),4,0);
        capsHbO2d=caps2D(HbO(end,find(lev==4)),4,0);
        capsHbO2d0=caps2D(HbO(1,find(lev==4)),4,0);
        capsHb2d=caps2D(Hb(end,find(lev==4)),4,0);
        capsHb2d0=caps2D(Hb(1,find(lev==4)),4,0);
        %% Calculating and Plotting fractional change
        delHbT= 100*(capsHbT2d-capsHbT2d0)./capsHbT2d0;
        delHbO= 100*(capsHbO2d-capsHbO2d0)./capsHbO2d0;
        delHb= 100*(capsHb2d-capsHb2d0)./capsHb2d0;
        index=find(lev==4);
        figure(1);
        ax(loopnr)= subplot(6,5,loopnr);
        delHbT=round(delHbT,2);
        delHbO=round(delHbO,2);
        delHb=round(delHb,2);
        imagesc(delHbT);axis off; 
        h=colorbar;caxis([-0.1 0.2]);
        set(h, 'Position', [.9 .11 .02 .8150]); %[left, bottom, width, height]
        pos=get(ax(loopnr), 'Position');
        set(ax(loopnr), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
        title([num2str(g/1000),'s']);
       
        figure(2);  ax(loopnr)=subplot(6,5,loopnr);
        imagesc(delHbO); axis off;
         h=colorbar;caxis([-0.1 0.2]);
        set(h, 'Position', [.9 .11 .02 .8150]); %[left, bottom, width, height]
        pos=get(ax(loopnr), 'Position');
        set(ax(loopnr), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
        title([num2str(g/1000),'s']);
        
        figure(3);  ax(loopnr)=subplot(6,5,loopnr);
        imagesc(delHb);axis off;
         h=colorbar;caxis([-0.1 0.15]);
        set(h, 'Position', [.9 .11 .02 .8150]); %[left, bottom, width, height]
        pos=get(ax(loopnr), 'Position');
        set(ax(loopnr), 'Position', [pos(1) pos(2) 0.85*pos(3) pos(4)]);
        title([num2str(g/1000),'s']);
        g=g+200;
        %   [delF,delV,delS,delHbT,delHbO,delHb]=percnt_chng(Fin1,V1,S1,HbT,HbO,Hb,X,Y,lev,bndry,j,CMRO2n,BETA,CMRO2n0); pause(0.1);
        
    end
end
figure(1);suptitle('HbT');
figure(2);suptitle('HbO');
figure(3);suptitle('Hb');

