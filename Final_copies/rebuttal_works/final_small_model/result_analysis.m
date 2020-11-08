close all;
clc;
clear all;
% load('output_4_plot_c3_amp=100.mat');
% disp('load the saved mat file');
% keyboard; 
inps=[2 10 80 90 200 250 350 400 600 1000];
for num=1:10
% filename = sprintf('r%d.mat',inps(num)) ;
filename = 'output_4_plot';
load(filename);
k=size(x,1);
% if k>20
% nm=1:5:k;
% else
     nm=1:k;
% end
load('XY.mat');
load('vessel_charac.mat'); 
load('lissomtrain1.mat');
% load('Z_constraint.mat');
file_neural='lissomtrain1.mat'; 
% [neu_size]=(size(Z_constraint));
[lissom,opts] = define(neu_size);

[lissom] = initialise4usage(lissom,file_neural); 
input1=input;
neu_length=neu_size(1)*neu_size(2);
eA=bndry(5); eV=bndry(6);



 nA=sum(layers(1:3));
nC=layers(4);
nV=sum(layers(5:7));
beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nA+nC);
CMRO2n0=x(1,3*nE+neu_length+nA+nC+1:end);
caps=sqrt(nC);
area_select=zeros(caps,caps);

area_select(2:4,2:4)=1;
neu_select=zeros(neu_size);
neu_select(2:4,2:4)=1;
lp=1;
ATPref=2.5;
In=1:k;
G=0:0.1:5;
Gt=3.8;
i=find(G==Gt);
g=Gt;
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
% Fin1=zeros(2,size(Fin,1));
% Fout1=zeros(2,size(Fin,1));
% Fin1(1,:)=Fin0'; Fout1(1,:)=Fout0';
% Fin1(2,:)=Fin';Fout1(2,:)=Fout';
% F=(Fin1+Fout1)/2;
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
% if (g<1) || (g>2)
%     input=zeros(size(input1));
% else
     input=input1;
% end
[neu,lissom]=inst_response_updat_depending_vsclr_fdk(lissom,opts,input,ATP,0);
% capsHb2d=caps2D(HbO(end,find(lev==4)),4,0);
%   capsHb2d0=caps2D(HbO(1,find(lev==4)),4,0);
% %% Calculating and Plotting fractional change
%   figure(j);subplot(331);imagesc(input); title(num2str(g));
%     meancmro2(lp)=mean(CMRO2n);
%         meanatp(lp)=mean(ATP);
%         meanF(lp)=mean(F(2,eA+1:eV-1));
%         meanpo2(lp)=mean(PO2tissn);
% %         meanoe(lp)=mean(OE1);
% %        meanpo2vess(lp)= mean(PO2_new);
%        meansat(lp)=mean(S(eA+1:eV-1));
% %        meanN(lp)=mean(reshape(neu,1,64));
%        meanbeta(lp)=mean(beta(eA+1:eV-1));
%        MHbT(lp)=mean(HbT(2,eA+1:eV-1)./HbT(1,find(lev==4)));
%        MHbO(lp)=mean(HbO(2,eA+1:eV-1)./HbO(1,find(lev==4)));
%        MHb(lp)=mean(Hb(2,eA+1:eV-1)./Hb(1,find(lev==4)));
%        lp=lp+1;
% end
delHbT= 100*(HbT-repmat(HbT(1,:),size(V1,1),1))./repmat(HbT(1,:),size(V1,1),1);
delHbO= 100*(HbO-repmat(HbO(1,:),size(V1,1),1))./repmat(HbO(1,:),size(V1,1),1);
delHb= 100*(Hb-repmat(Hb(1,:),size(V1,1),1))./repmat(Hb(1,:),size(V1,1),1);
capsHb2d=caps2D(delHbT(end,find(lev==4)),4,0);
caps_final=capsHb2d.*area_select./max(max(capsHb2d.*area_select));
 Avg_neu(num)=sum(sum(neu.*neu_select))/numel(find(neu_select>0));
Avg_HbT(num)=sum(sum(caps_final))/numel(find(area_select>0));

end
plot(Avg_neu); hold on ; plot(Avg_HbT,'r');