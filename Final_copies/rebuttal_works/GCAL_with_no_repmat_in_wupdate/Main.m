clear all;close all;clc;
load X_test_4vasc_amp5.mat;
tic;
a=X1 ;
% a=a1(:,:,1:3);
xdim = size(a(:,:,1));
lissom_size=[32,32];

[neu_size]=lissom_size;
%% intialization
[lissom,opts] = define(neu_size);                      %fn1
proceed=check_sheet_size_rf(lissom,a);
lissom = initialise_new(lissom,opts,X1);  
ATPmax=10;%(opts.p+op4ts.q-opts.r)/2;
ATPref=2.5;
ATP=ATPref*ones(lissom.layers{2}.dim(1)*lissom.layers{2}.dim(2),1);%fn3
%% training
ijk=1;
Conv=0;
tic;
if proceed(1,1)==1 && proceed(1,2)==1
    for ww=1:opts.epoch                     %30 epochs to train the 90 inputs per epoch randomly selected
%         fprintf('\nepoch %d ',ww); %ww-epoch i-input
        p=randperm(size(a,3));              %90 random orientations stored in p input to train
        for i= 1: size(a,3)                 %for each input orientation
             inpsel=a(:,:,p(i));            %the image for the input orientation is taken from a
            fprintf('epoch %d input %d\n',ww,i);
                for jj = 2:numel(lissom.layers)                  %for each layer, here 2 layers i/p and V1           %fn4
                    if jj==2
                        INPUT=reshape(inpsel,numel(inpsel),1);
%                          figure(11);subplot(121);imagesc(reshape(INPUT,sqrt(size(INPUT,1)),sqrt(size(INPUT,1))));title('Input');colorbar;pause(0.01);                        
                        [lissom]= RF_master_new(lissom,jj,INPUT);
                    end
%                     if jj==3
%                         INPUT=inpsel1;%reshape(inpsel1(:,:,j),numel(inpsel1(:,:,j)),1);
%                         [lissom]= RF_master_new(lissom,jj,INPUT);%select receptive field for each neuron --'X'       %fn5
%                     end
                    [lissom] = lat_dynamics_final(lissom,jj,opts,lissom.layers{jj}.padded_input,ATP,ATPmax); %lateral dynamics
%                     if jj==2
%                         inpsel1=lissom.layers{jj}.Zold;
%                     end
                    [lissom] = wupdate(lissom,jj,opts);%update weights            %fn8
                    
                    ijk=ijk+1;
                end
                lissom.layers{jj}.Zprev = lissom.layers{jj}.Zold ;
            
            for ii=2:numel(lissom.layers)
                lissom.layers{ii}.Zprev =zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
                lissom.layers{ii}.Zold =zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
                lissom.layers{ii}.Z=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
            end
        end
%         if (ww==1)||(mod(ww,2)==0)
%             save_name = sprintf('Epoch%d.mat',ww);
%             save(save_name,'lissom','opts','label');
%         end
    end
else
    error('adjust rf');
end

%% saving
MAPS(lissom,opts,X1,label,ATP,ATPmax);
save('lissom_10atpmax_1p_10q_8r.mat','lissom','opts','ATPmax','ATPref');
% % cd F:\Github_team\GCAL_model\RF_PF_variations_oct_7_2020\optimize_param
%  str =  'lesioned' ;%+ string(opts.q(2))+'_R_' + string(opts.r(2))+'.mat';%'PF_' + string(lissom.layers{2}.rf_fb(1))+'.mat';
%  save(str,'lissom','opts','label');


% %% testing
% fprintf('\n');
% yo=toc;
% if yo>60 && yo<=3600
%     yt=yo/60;
%     ds=['elapsed time before map = ',num2str(yt),'min'];
%     disp(ds);
% elseif yo>3600
%     yt=yo/(60*60);
%     ds=['elapsed time before map = ',num2str(yt),'hours'];
%     disp(ds);
% end
% Xtest=X1;
% [lissom]=make_map_newDIavg(lissom,opts,Xtest,label);    %fn11