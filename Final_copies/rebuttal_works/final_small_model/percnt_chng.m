function [delF,delV,delS,delHbT,delHbO,delHb]=percnt_chng(F,V,S,HbT,HbO,Hb,X,Y,lev,bndry,ind,PO2tissn,beta,CMRO2n0)
eA=bndry(5); eV=bndry(6);
caps=sqrt(eV-eA-1);
delF= 100*(F-repmat(F(1,:),size(F,1),1))./repmat(F(1,:),size(F,1),1);
delV= 100*(V-repmat(V(1,:),size(V,1),1))./repmat(V(1,:),size(V,1),1);
delS= 100*(S-repmat(S(1,:),size(S,1),1))./repmat(S(1,:),size(S,1),1);
delHbT= 100*(HbT-repmat(HbT(1,:),size(V,1),1))./repmat(HbT(1,:),size(V,1),1);
delHbO= 100*(HbO-repmat(HbO(1,:),size(V,1),1))./repmat(HbO(1,:),size(V,1),1);
delHb= 100*(Hb-repmat(Hb(1,:),size(V,1),1))./repmat(Hb(1,:),size(V,1),1);
po2ini=CMRO2n0;
delPO2=100*((PO2tissn-po2ini)./po2ini);
delbeta=100*(beta-repmat(beta(1,:),size(beta,1),1))./repmat(beta(1,:),size(beta,1),1);
% 
% plot_func(delF(end,:)',X,Y);
% plot_func(delV(end,:)',X,Y);
% plot_func(delS(end,:)',X,Y);
% plot_func(delHbT(end,:)',X,Y);
% plot_func(delHbO(end,:)',X,Y);
% plot_func(delHb(end,:)',X,Y);

% plot_func(delF(end,eA+1:eV-1)',X(eA+1:eV-1),Y(eA+1:eV-1),'C');
% plot_func(delV(end,eA+1:eV-1)',X(eA+1:eV-1),Y(eA+1:eV-1),'C');
% plot_func(delS(end,eA+1:eV-1)',X(eA+1:eV-1),Y(eA+1:eV-1),'C');
% plot_func(delHbT(end,eA+1:eV-1)',X(eA+1:eV-1),Y(eA+1:eV-1),'C');
% plot_func(delHbO(end,eA+1:eV-1)',X(eA+1:eV-1),Y(eA+1:eV-1),'C');
% plot_func(delHb(end,eA+1:eV-1)',X(eA+1:eV-1),Y(eA+1:eV-1),'C');
% for ii=1:size(delV,2)
%     if delV(2,ii)<1
%         delV(2,ii)=-10;
%     end
%     if delHbT(2,ii)<1
%        delHbT(2,ii)=-10;
%     end
%  if delHb(2,ii)<10
%        delHb(2,ii)=-10;
%     end
% end

for i=4%2:lev(end)-1

    j=find(lev==i);
%       figure(ind); subplot(3,3,2)
%       plotcaps(delbeta(end,:),i,0);title('BETA')

   figure(ind)
    
   subplot(3,3,3); 
   plotcaps(delF(end,j),i,0); title('F')
       figure(ind);
     subplot(3,3,4);
    plotcaps(delV(end,j),i,0);title('V') % change this value 8 according to no. of cappilliaries. now 64 caps, so sqrt 8 taken
 
       figure(ind);
     subplot(3,3,5);
    plotcaps(delHbT(end,j),i,0); title('HbT')
   figure(ind); subplot(3,3,6);
    plotcaps(delHbO(end,j),i,0); title('HbO')
   figure(ind); subplot(3,3,7);
    
    plotcaps(delHb(end,j),i,0); title('Hb')
   figure(ind); subplot(3,3,8);
    plotcaps(delS(end,j),i,0); title('S')
    
%      figure(ind); subplot(3,3,9);
%     plotneuside(delPO2,8,8); title('CMRO2n')
end