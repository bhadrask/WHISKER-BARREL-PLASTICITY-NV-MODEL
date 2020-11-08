function [lissom,label]= MAPS(lissom,opts,X1,label,ATP,ATPmax) 
xx=[];
label(1,9)=-1;
opts.niter=1;
for i=1:size(X1,3)     %testing for 8 input videos
    for ijk=2:numel(lissom.layers)
        lissom.layers{ijk}.ZL2 =zeros(lissom.layers{ijk}.dim(1)*lissom.layers{ijk}.dim(2),1);
        lissom.layers{ijk}.Zold =zeros(lissom.layers{ijk}.dim(1)*lissom.layers{ijk}.dim(2),1);
    end
    inpsel=(X1(:,:,i));  %taking 1 input image at a time inpsel=cell2mat(a(:,p(i)));
    steps=size(inpsel,3);
    for j=1:steps
        for jj=2:numel(lissom.layers)
                                       figure(1);subplot(121);imagesc(inpsel(:,:,j));colorbar;pause(0.1);
            if jj==2
                INPUT=reshape(inpsel(:,:,j),numel(inpsel(:,:,j)),1);
            end
            if jj==3
                INPUT=reshape(inpsel2(:,:,j),numel(inpsel2(:,:,j)),1);
            end
            [lissom]= RF_master_new(lissom,jj,INPUT);
            
            [lissom] = lat_dynamics_final_4_response(lissom,jj,opts,lissom.layers{jj}.padded_input,ATP,ATPmax); %lateral dynamics
            if jj==2
                inpsel2(:,:,j)=lissom.layers{jj}.Zold;
                zim=reshape(lissom.layers{jj}.Zold,lissom.layers{jj}.dim);
                                          figure(1);        subplot(122);imagesc(zim);colorbar;pause(0.1);
            end
            if jj==3
                zim=reshape(lissom.layers{jj}.Zold,lissom.layers{jj}.dim);
%                                                  subplot(133);imagesc(zim);colorbar;pause(0.01);
            end
            if j==size(inpsel,3)
                Z2=reshape(lissom.layers{jj}.Zold,lissom.layers{jj}.dim);
                xx=cat(3,xx,Z2);
            end
            lissom.layers{jj}.converg(i,j)=(sum(sum(lissom.layers{jj}.ZL2-lissom.layers{jj}.Zold).^2))/lissom.layers{2,1}.dim(1,1)/lissom.layers{2,1}.dim(1,1);
            lissom.layers{jj}.ZL2=lissom.layers{jj}.Zold;
        end
    end
end
[gg,q]=max(xx,[],3);
% figure; imagesc(labeldir(q));title('layer2 direction label');colorbar;
loc= gg<=0;
q(loc)=length(label);

     Dwhisker=numel(find(label(q)<=30 ));
     number0=numel(find(label(q)<=0));
      Dwhisker= Dwhisker-number0
Cwhisker=numel(find(label(q)>30))
fraction=Dwhisker/Cwhisker
figure(5); imagesc(label(q));title('barrels');colorbar;