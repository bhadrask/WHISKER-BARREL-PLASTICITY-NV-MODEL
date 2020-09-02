function [yy,lissom]=inst_response_updat_depending_vsclr_fdk(lissom,opts,X1,ATP,tr,ATPmax,thr)
yy=[];
ATP1=reshape(ATP,lissom.layers{2}.dim(1), lissom.layers{2}.dim(2));
inpsel=X1;

for jj=2:numel(lissom.layers)
    [lissom] = rf_m(lissom,double(inpsel),jj,opts);
    [lissom] = activate(lissom,jj,opts);
    
    [lissom] = lat_dynamics_final_vas_intr(lissom,jj,opts,ATP1,ATPmax,thr);
    
    yy=lissom.layers{jj}.Zold;
    
    
end
% figure(2);
% subplot(211);imagesc(inpsel);subplot(212);imagesc(yy); pause(0.001);
if tr==1
for kk=2:numel(lissom.layers)
            [lissom] = wupdate(lissom,kk,opts);%update weights
            [lissom]=make_unit_norm(lissom,kk,1);
            lissom.layers{kk}.Zold = zeros(lissom.layers{kk}.dim);
            lissom.layers{kk}.Z = zeros(lissom.layers{kk}.dim);% reset activation to zero before new input
end
end