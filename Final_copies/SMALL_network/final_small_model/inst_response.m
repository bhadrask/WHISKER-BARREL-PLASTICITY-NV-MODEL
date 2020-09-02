function [yy,MRES]=inst_response(lissom,opts,X1,label,ATP,ATPmax,thr)
yy=[];xx=[];xx2=[];MRES=0;xx3=[];
label(length(label)+1)=-1;
ATP1=reshape(ATP,lissom.layers{2}.dim(1), lissom.layers{2}.dim(2));
for i=1:size(X1,3)
    inpsel=X1(:,:,i);
    figure(i); subplot(2,1,1);imagesc(inpsel);
    for jj=2:numel(lissom.layers)
        [lissom] = rf_m(lissom,double(inpsel),jj,opts);
        [lissom] = activate(lissom,jj,opts);
        
        %     [lissom] = lat_dynamics(lissom,jj,opts,Z_constraint);
           [lissom] = lat_dynamics_no_sigmoid(lissom,jj,opts,ATP1,ATPmax,thr);
%         [lissom]=lat_dynamics_final_vas_intr(lissom,jj,opts,ATP1,ATPmax);
        
        yy=lissom.layers{jj}.Zold;
        yy1=lissom.layers{jj}.Z;
           yy2=lissom.layers{jj}.Z_nosig;
        xx=cat(3,xx,lissom.layers{jj}.Zold);
        xx2=cat(3,xx2,lissom.layers{jj}.Z);
               xx3=cat(3,xx3,lissom.layers{jj}.Z_nosig);
        for kk=2:numel(lissom.layers)
            lissom.layers{kk}.Zold = zeros(lissom.layers{kk}.dim);
            lissom.layers{kk}.Z = zeros(lissom.layers{kk}.dim);% reset activation to zero before new input
        end
    end
    meanresp=max(max(yy));
    MRES=MRES+meanresp;
    figure(i);
    subplot(2,1,2);imagesc(yy1);title(['input',num2str(i)]); pause(0.5)
    
end;
MRES=MRES/size(X1,3);
%   [gg,q]=max(xx,[],3)
%    % figure; imagesc(labeldir(q));title('layer2 direction label');colorbar;
%        loc=find(gg<=0);
%      q(loc)=length(label);
%     figure(); imagesc(label(q));title('barrels');colorbar;
for kk=1:size(xx,3)
    lr=label(kk)*ones(size(xx,1),size(xx,1));
    mp(:,:,kk)=lr.*xx(:,:,kk);
end
sigmap=round(sum(mp,3)./sum(xx,3));
sigmap(isnan(sigmap))=-10;

[gg,q]=max(xx,[],3);
%   MLFP=mean(LFP);
% figure; imagesc(labeldir(q));title('layer2 direction label');colorbar;
loc=find(gg<=0);
q(loc)=length(label);
figure(); imagesc(label(q));title('barrels');colorbar;
[gg3,q3]=max(xx3,[],3);
Dwhisker=numel(find(label(q3)<=30 ));
number0=numel(find(label(q3)<=0));
Dwhisker= Dwhisker-number0
Cwhisker=numel(find(label(q3)>30))
fraction=Dwhisker/Cwhisker
figure(); imagesc(label(q3));title(['barrels no sigmoid, D/C ratio=',num2str(fraction),', Mean Activity=',num2str(MRES)]);colorbar;
%      figure();imagesc(sigmap);title('weighted average map');

end

