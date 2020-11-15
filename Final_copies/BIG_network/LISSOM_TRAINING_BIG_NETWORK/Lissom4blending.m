clear all;close all;clc;
% load Final_X_test_sigma_6.mat;
load Final_X_test_sigma_4.mat
tic;
a=X1 ;xdim = size(a(:,:,1));
%% section for lissom constraining
lissom_size=[64,64];
% A=3;
% [Z_constraint] = get_constraint(A,lissom_size);
load('Z_constraint');
[neu_size]=(size(Z_constraint));
%% 
[lissom,opts] = define(neu_size);
 lissom = initialise(lissom,xdim,a,opts);
% lissom = initialise_old(lissom,xdim,a,opts);
ATPmax=10;%(opts.p+opts.q-opts.r)/2;
ATP=zeros(lissom.layers{2}.dim(1), lissom.layers{2}.dim(2));
for ini=2: numel(lissom.layers)
[lissom]=make_unit_norm(lissom,ini,1);% the last factor decides the norm of normalization , either 1or 2
end
for ww=1:opts.epoch
  
    p=randperm(size(a,3));
    for i= 1: size(p,2) %for each input
         
        inpsel=a(:,:,p(i));
        fprintf('epoch %d input %d\n',ww,i);
%         figure(1) ;%subplot(3,3,1);
%         imagesc(inpsel);title(num2str(max(max(inpsel))))
        for jj = 2 : numel(lissom.layers) %for each layer
            
            [lissom] = rf_m(lissom,inpsel,jj,opts);%select receptive field for each neuron --'X'
            [lissom] = activate(lissom,jj,opts);% initial activity from previous layer input-- 'Z'
            
            if jj>1 && jj<numel(lissom.layers) && opts.control==1
                [lissom] = feedback(lissom,jj,opts);
            end
             [lissom] = lat_dynamics_final(lissom,jj,opts,Z_constraint,ATP,ATPmax); %lateral dynamics
%              [lissom] =lat_dynamics_old(lissom,jj,opts,Z_constraint);
            inpsel=lissom.layers{jj}.Zold;
            
        end
        for kk=2:numel(lissom.layers)
            [lissom] = wupdate(lissom,kk,opts);%update weights
            [lissom]=make_unit_norm(lissom,kk,1);
            lissom.layers{kk}.Zold = zeros(lissom.layers{kk}.dim);
            lissom.layers{kk}.Z = zeros(lissom.layers{kk}.dim);% reset activation to zero before new input
        end
    end
end
lissomaff= lissom.layers{2}.waff;
lissomexc= lissom.layers{2}.wnbrexc;
lissominhb= lissom.layers{2}.wnbrinhb;
if opts.control==1
    lissomfb=lissom.layers{2}.wfdbk;
    
else
    lissomfb=0;
end
%rf1=lissom.layers{2}.rf;

save('lissom_10atpmax_newopts_14q_8r_sigma4.mat','lissom','opts','ATPmax');
% gama1=opts.p;gama2=opts.q;gama3=opts.r;dim1=lissom.layers{2}.dim;rf1=lissom.layers{2}.rf;
% save('lissomtrain1_10atpmax.mat','lissomaff','lissomexc','lissominhb','lissomfb','gama1','gama2','gama3','dim1','rf1');
% load('Final_X_test_sigma_6_AMP5.mat');
Xtest=X1;
inst_response(lissom,opts,Xtest,Z_constraint,label,ATP,ATPmax);
yo=toc;
if yo>60 && yo<=3600
     yt=yo/60;
     ds=['elapsed time = ',num2str(yt),'min'];
    disp(ds);
elseif yo>3600
    yt=yo/(60*60);
    ds=['elapsed time = ',num2str(yt),'hours'];
     disp(ds);
end
    