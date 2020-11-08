% clear all;
close all;
clc;
yy=[];
xx=[];
neu_length=neu_size(1)*neu_size(2);
for i=1:size(X1,3)
    input=X1(:,:,i);
    ATP=2.5*(ones(1,neu_length));%x(end,3*nE+1:3*nE+neu_length);
[neu,lissom]=inst_response_updat_depending_vsclr_fdk(lissom,opts,input,ATP,0);
figure();imagesc(neu);pause(0.5);
   yy=neu;
df=find(neu>0.1);
numcell(i)=numel(df);
if i==1 || i==3
    neu(neu>0.46)=1.5;
end
       xx=cat(3,xx,neu);
       
end
  [gg,q]=max(xx,[],3);
    figure(); imagesc(label(q));title('barrels');colorbar;
  for kk=1:size(xx,3)
      lr=label(kk)*ones(size(xx,1),size(xx,1));
      mp(:,:,kk)=lr.*xx(:,:,kk);
  end
sigmap=(sum(mp,3)./sum(xx,3));
sigmap(isnan(sigmap))=-10;
    figure();imagesc(sigmap);title('weighted average map');
    ratio=sum(numcell(1:3))/sum(numcell(4:6))