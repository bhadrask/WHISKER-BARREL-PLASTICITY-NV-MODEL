clear all;close all;clc;
load X_test.mat;
tic;
a=X1 ;xdim = size(a(:,:,1));
%% section for lissom constraining
lissom_size=[8,8];

[neu_size]=lissom_size;
%% 
[lissom,opts] = define(neu_size);
 lissom = initialise(lissom,xdim,a,opts);
% lissom = initialise_old(lissom,xdim,a,opts);
ATPmax=10;%(opts.p+op4ts.q-opts.r)/2;
ATPref=2.5;
ATP=ATPref*ones(lissom.layers{2}.dim(1), lissom.layers{2}.dim(2));
figure(3);imagesc(cell2mat(lissom.layers{2}.waff));
for ini=2: numel(lissom.layers)
[lissom]=make_unit_norm(lissom,ini,1);% the last factor decides the norm of normalization , either 1or 2
end
for ww=1:opts.epoch
  
    p=randperm(size(a,3));
    for i= 1: size(p,2) %for each input
         
        inpsel=a(:,:,p(i));
        fprintf('epoch %d input %d\n',ww,i);
        figure(1) ;%subplot(3,3,1);
        imagesc(inpsel);title(num2str(max(max(inpsel))))
        for jj = 2 : numel(lissom.layers) %for each layer
            
            [lissom] = rf_m(lissom,inpsel,jj,opts);%select receptive field for each neuron --'X'
            [lissom] = activate(lissom,jj,opts);% initial activity from previous layer input-- 'Z'
            
           

            [lissom] = lat_dynamics_final(lissom,jj,opts,ATP,ATPmax); %lateral dynamics

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


save('lissom_10atpmax_1p_10q_8r.mat','lissom','opts','ATPmax','ATPref');

Xtest=X1;
inst_response(lissom,opts,Xtest,label,ATP,ATPmax);
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
    
figure(4);imagesc(cell2mat(lissom.layers{2}.waff));