function [yy,lissom]=inst_response_updat_depending_vsclr_fdk(lissom,opts,X1,ATP,tr,ATPmax,thr)
yy=[];
ATP1=ATP;
inpsel=X1;
if numel(lissom.layers)>2
    error('not coded for 3 layer');
end
jj=2;
INPUT=reshape(inpsel,numel(inpsel),1);
[lissom]= RF_master_new(lissom,jj,INPUT);
[lissom] = lat_dynamics_final_vas_intr(lissom,jj,opts,lissom.layers{jj}.padded_input,ATP1,ATPmax,thr); %lateral dynamics
%
yy=lissom.layers{jj}.Zold;

% figure(2);
% subplot(211);imagesc(inpsel);subplot(212);imagesc(yy); pause(0.001);
if tr==1
    ii=2;
    [lissom] = wupdate(lissom,ii,opts);%update weights
    lissom.layers{ii}.Zprev =zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
    lissom.layers{ii}.Zold =zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
    lissom.layers{ii}.Z=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
    
end
end