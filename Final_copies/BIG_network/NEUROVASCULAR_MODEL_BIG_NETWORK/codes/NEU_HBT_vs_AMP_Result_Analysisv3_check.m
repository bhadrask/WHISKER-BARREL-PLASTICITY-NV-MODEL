% close all;
clc;
clear all;

for num=1:1;
load('Final_amp200.mat');
    load('XY.mat');
    load('vessel_charac.mat');
    input1=IN.input;
    
    neu_length=neu_size(1)*neu_size(2);
    eA=bndry(5); eV=bndry(6);
    
    g=-0.1;
    nC=layers(4);
    
    beta0=x(1,3*nE+neu_length+1:3*nE+neu_length+nE);
    CMRO2n0=x(1,3*nE+neu_length+nE+1:end);
    V0=x(1,1:nE);
    S0=x(1,nE+1:2*nE);
    
    caps=sqrt(nC);
%     area_select1=zeros(caps,caps);
%     area_select1(11:48,20:50)=1;
%     neu_select=zeros(neu_size);

%     neu_select(34:36,30:32)=1;
    
     area_select1=zeros(caps,caps);
% %     area_select1(25:37,29:42)=1;
%      area_select1(26:40,26:40)=1;
area_select1(30:36,32:38)=1;
    neu_select=zeros(neu_size);
    neu_select(29:33,33:37)=1;
      
    
    
    lp=1;
    ATPref=2.5;
    
    G=0:0.1:5;
    Gt=2.5;
    j=find(G==Gt);
    g=Gt;
    for i=j
        
        
        V=x(i,1:nE);
        S=x(i,nE+1:2*nE);
        
        beta=x(i,3*nE+neu_length+1:3*nE+neu_length+nE);
        
        V1=zeros(2,size(V0,2));
        S1=zeros(2,size(V0,2));
        
        V1(1,:)=V0;
        V1(2,:)=V';
        S1(1,:)=S0;
        S1(2,:)=S';
        HbT=2.3.*S1.*V1;%2.3*V1;
        HbO=2.3.*S1.*V1;
        Hb=2.3*(1-S1).*V1;
        
        input=input1;
        [~,idx]=min(abs(Ttime-Gt+0.6));
        gtime=Ttime(idx);
        neu=neuout(:,:,idx);
        capsHb2d=caps2D(HbT(end,lev==4),4,0);
        capsHb2d0=caps2D(HbT(1,lev==4),4,0);
       Beta0=caps2D(beta0(lev==4),4,0);
       Beta=caps2D(beta(lev==4),4,0);
        HbT_var_2d_1=(capsHb2d-capsHb2d0)./capsHb2d0;
        HbT_var_2d_1(abs(HbT_var_2d_1)<1e-4)=0;
        if max(max(HbT_var_2d_1))==0
            HbT_var_2d=0;
        else
        HbT_var_2d= area_select1.*HbT_var_2d_1/max(max(HbT_var_2d_1));
        end
        Avg_Hbt_var=sum(sum(HbT_var_2d)/numel(find(area_select1>0)));
        sum(sum(HbT_var_2d))
        AA=(capsHb2d-capsHb2d0)./capsHb2d0;
        percentage=max(max(AA))
        figure();subplot(211);imagesc((capsHb2d-capsHb2d0)./capsHb2d0); title('HbO')
        subplot(212);imagesc(neu); 
        figure();imagesc((Beta-Beta0)./Beta0);title('beta')
        figure;surf((Beta-Beta0)./Beta0)
        if max(max(neu.*neu_select))==0
            neu1=0;
        else
            neu1=neu.*neu_select;%./max(max(neu.*neu_select));
        end
      
        Avg_neu(num)=sum(sum(neu1))/numel(find(neu_select>0));
        Avg_HbT(num)=Avg_Hbt_var;
        HHbt=caps2D(HbT(end,lev==4),4,0);
        %         Hmax(num)=max(max(delHbT));
    end
end

% figure;plot(Avg_neu); hold on ; plot(Avg_HbT,'r');

