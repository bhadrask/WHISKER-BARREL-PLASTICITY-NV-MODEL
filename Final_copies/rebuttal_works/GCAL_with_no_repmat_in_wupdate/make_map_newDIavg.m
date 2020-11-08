function [lissom]=make_map_newDIavg(lissom,opts,X1,label)
lissom.layers{2, 1}.xx=[];lissom.layers{2, 1}.xx8=[];
lissom.layers{3, 1}.xx=[];lissom.layers{3, 1}.xx8=[];
label(1,9)=-1;
opts.niter=1;
for i=1:size(X1,2)     %testing for 8 input videos
    for ijk=2:numel(lissom.layers)
        lissom.layers{ijk}.ZL2 =zeros(lissom.layers{ijk}.dim(1)*lissom.layers{ijk}.dim(2),1);
        lissom.layers{ijk}.Zold =zeros(lissom.layers{ijk}.dim(1)*lissom.layers{ijk}.dim(2),1);
    end
    inpsel=cell2mat(X1(:,i));  %taking 1 input image at a time inpsel=cell2mat(a(:,p(i)));
    steps=size(inpsel,3);
    for j=1:steps
        for jj=2:numel(lissom.layers)
%                                       figure(12);subplot(131);imagesc(inpsel(:,:,j));colorbar;pause(0.01);
            if jj==2
                INPUT=reshape(inpsel(:,:,j),numel(inpsel(:,:,j)),1);
            end
            if jj==3
                INPUT=reshape(inpsel2(:,:,j),numel(inpsel2(:,:,j)),1);
            end
            [lissom]= RF_master_new(lissom,jj,INPUT);
            
            [lissom] = lat_dynamics_final(lissom,jj,opts,lissom.layers{jj}.padded_input); %lateral dynamics
            if jj==2
                inpsel2(:,:,j)=lissom.layers{jj}.Zold;
                zim=reshape(lissom.layers{jj}.Zold,lissom.layers{jj}.dim);
%                                                   subplot(132);imagesc(zim);colorbar;pause(0.01);
            end
            if jj==3
                zim=reshape(lissom.layers{jj}.Zold,lissom.layers{jj}.dim);
%                                                  subplot(133);imagesc(zim);colorbar;pause(0.01);
            end
            if j==size(inpsel,3)
                Z2=reshape(lissom.layers{jj}.Zold,lissom.layers{jj}.dim);
                lissom.layers{jj}.xx=cat(3,lissom.layers{jj}.xx,Z2);
            end
            lissom.layers{jj}.converg(i,j)=(sum(sum(lissom.layers{jj}.ZL2-lissom.layers{jj}.Zold).^2))/lissom.layers{2,1}.dim(1,1)/lissom.layers{2,1}.dim(1,1);
            lissom.layers{jj}.ZL2=lissom.layers{jj}.Zold;
        end
    end
end
for some = 2:numel(lissom.layers)
    clear xxnorm
    clear mp_cos
    clear mp_sin
    clear add_cos
    clear add_sin
    clear add_xxnorm
    xxnorm=(lissom.layers{some}.xx)./repmat(max(lissom.layers{some}.xx,[],3),[1,1,size(lissom.layers{some}.xx,3)]);
    xxnorm(isnan(xxnorm)) = 0;
       
    [gg,q]=max(xxnorm,[],3);%q = layer which has the max value
    %     [ii jj] = find(isnan(gg));
    %     for somenum=1:size(ii,1)
    %         gg(ii(somenum,1),jj(somenum,1))=-20;
    %     end
    loc=find(gg<=0);
    q(loc)=length(label);
    
    if any(find(q==9))
        colors = [0 0 0; % black for -1
        1 0 0; % red for 0
        1 0.5 0; % orange for 1
        1 1 0; % yellow for 2
        0 1 0; % green for 3, etc.
        0.5 0 1; %purple for 4
        0 0.2 1; % blue for 5
        0.2 0.8 1; % cyan for 6
        1 0 1];% pink for 7
    else
        colors = [%0 0 0; % black for -1
        1 0 0; % red for 0
        1 0.5 0; % orange for 1
        1 1 0; % yellow for 2
        0 1 0; % green for 3, etc.
        0.5 0 1; %purple for 4
        0 0.2 1; % blue for 5
        0.2 0.8 1; % cyan for 6
        1 0 1];% pink for 7
    end
    
    %the index which has the negative value in gg, for that index in q = 5
    figure();%subplot(121);
    if some==2
    imagesc(label(q));title(['Barrels of Neural layer']);colorbar;% 1 2 3 4 5 in q is replaced by 180 135 90 45 -1 (label) respectively
    end
    if some==3
    imagesc(label(q));title(['Barrels of Vascular layer ']);colorbar;% 1 2 3 4 5 in q is replaced by 180 135 90 45 -1 (label) respectively
    end
    colormap(colors);
    % subplot(122);histogram(label(q));title('distribution layer 1');
    
    
    for kk=1:4%size(lissom.layers{2}.xx8,3)
        lr=2*label(kk)*ones(lissom.layers{some, 1}.dim(1,1),lissom.layers{some, 1}.dim(1,1)); % 10 x 10 full of 180 135 90 45
        mp_cos(:,:,kk)=cosd(lr).*xxnorm(:,:,kk);                % 10 x 10 x 1 mp=lr*xx
        mp_sin(:,:,kk)=sind(lr).*xxnorm(:,:,kk);
    end
    for pi=1:lissom.layers{some, 1}.dim(1,1)
        for pj=1:lissom.layers{some, 1}.dim(1,1)
            add_cos(pi,pj)=sum(mp_cos(pi,pj,:));
            add_sin(pi,pj)=sum(mp_sin(pi,pj,:));
            lissom.layers{some}.sigmap(pi,pj)=2*(atand(add_cos(pi,pj)./(add_sin(pi,pj))));  % orientation preference OP
            Sety(pi,pj)=sqrt(add_cos(pi,pj).^2 + add_sin(pi,pj).^2);
            add_xxnorm(pi,pj)=sum(xxnorm(pi,pj,1:4));
            lissom.layers{some}.Selectivity(pi,pj)=Sety(pi,pj)./add_xxnorm(pi,pj);    % orientation selectivity
            if lissom.layers{some}.sigmap(pi,pj)<=0
                lissom.layers{some}.sigmap(pi,pj)=180 + lissom.layers{some}.sigmap(pi,pj);              % orientation preference
            end
        end
    end
    
    [gg,q]=max(xxnorm,[],3);%q = layer which has the max value
    for num_x=1:lissom.layers{some, 1}.dim(1,1)
        for num_y=1:lissom.layers{some, 1}.dim(1,1)
            if q(num_x,num_y)<=4
                q_opp(num_x,num_y)=q(num_x,num_y)+4;
                gg_opp(num_x,num_y)=xxnorm(num_x,num_y,q_opp(num_x,num_y));
                lissom.layers{some}.di(num_x,num_y)=1-gg_opp(num_x,num_y)/gg(num_x,num_y);
            end
            if q(num_x,num_y)<=8 && q(num_x,num_y)>4
                q_opp(num_x,num_y)=q(num_x,num_y)-4;
                gg_opp(num_x,num_y)=xxnorm(num_x,num_y,q_opp(num_x,num_y));
                lissom.layers{some}.di(num_x,num_y)=1-gg_opp(num_x,num_y)/gg(num_x,num_y);
            end
        end
    end
    
    rf = lissom.layers{2}.rf_fb(1,1); st = lissom.layers{3}.stride;
    groupn = (lissom.layers{some, 1}.dim(1,1) - rf)/st + 1;
    for somenum = 1
        p=1;
        for l = 1 :groupn%lissom.layers{some}.dim(1)     %1 to 10 %[4:9] to [1:6], [4:9] to [4:9], [4:9] to [7:12]....
            m=1;
            for k = 1 :groupn%lissom.layers{some}.dim(2) %1 to 10 %[1:6] to [1:6], [4:9] to [1:6], [7:12] to [1:6]....            /max(xx(l,k,:));
                lissom.layers{some,1}.selmat1{l,k} = lissom.layers{some, 1}.sigmap(p:p+(rf-1),m:m+(rf-1),somenum);
                lissom.layers{some,1}.selmat2{l,k} = lissom.layers{some, 1}.Selectivity(p:p+(rf-1),m:m+(rf-1),somenum);
                lissom.layers{some,1}.selmat3{l,k} = lissom.layers{some, 1}.di(p:p+(rf-1),m:m+(rf-1),somenum);
                if  some == 2 % average of neural activity
                    lissom.layers{some,1}.sigmapavg(l,k,somenum) = sum(sum((lissom.layers{some,1}.selmat1{l,k})))/ (rf*rf);
                    lissom.layers{some,1}.Selectivityavg(l,k,somenum) = sum(sum((lissom.layers{some,1}.selmat2{l,k})))/ (rf*rf);
                    lissom.layers{some,1}.diavg(l,k,somenum) = sum(sum((lissom.layers{some,1}.selmat3{l,k})))/ (rf*rf);                                        
                end
                if  some ==3  % sum of the vascular activity
                    lissom.layers{some,1}.sigmapavg(l,k,somenum) = lissom.layers{some, 1}.sigmap(l,k,somenum);
                    lissom.layers{some,1}.Selectivityavg(l,k,somenum) = lissom.layers{some, 1}.Selectivity(l,k,somenum);
                    lissom.layers{some,1}.diavg(l,k,somenum) = lissom.layers{some, 1}.di(l,k,somenum); 
                end
                m=k*st+1;   %q=1,4,7,11...
            end
            p=l*st+1;       %p=1,4,7,11...
        end
    end
    
end

figure();
sct=0;
for pi=1:groupn
    for pj=1:groupn
        diff=abs(lissom.layers{3}.sigmapavg(pi,pj)-lissom.layers{2}.sigmapavg(pi,pj));
        if (diff<60)
            scatter(lissom.layers{2}.sigmapavg(pi,pj),lissom.layers{3}.sigmapavg(pi,pj),'b','LineWidth',2);hold on;
            sct=sct+1;
            idx_pi(sct)=pi;
            idx_pj(sct)=pj;
        end
    end
end
plot(0:1:180,0:1:180,'r','LineWidth',2);
set(gca,'fontweight','bold','LineWidth',2)
xlabel('Orientation Preference of Neurons','fontweight','bold');ylabel('Orientation Preference of Blood Vessels','fontweight','bold');title('Orientation Preference','fontweight','bold');%ylim([0 1]);xlim([0 1]);

figure();
for pi=1:length(idx_pi)
    %          scatter(Selectivity(pi,pj),Selectivity1(pi,pj),'b');hold on;
    scatter(lissom.layers{2}.Selectivityavg(idx_pi(pi),idx_pj(pi)),lissom.layers{3}.Selectivityavg(idx_pi(pi),idx_pj(pi)),'b','LineWidth',2);hold on;
    sel_selectivity(pi,1)=lissom.layers{2}.Selectivityavg(idx_pi(pi),idx_pj(pi));
    sel_selectivity(pi,2)=lissom.layers{3}.Selectivityavg(idx_pi(pi),idx_pj(pi));
end
plot(0:0.1:1,0:.1:1,'r','LineWidth',2);
set(gca,'fontweight','bold','LineWidth',2)
xlabel('Orientation Selectivity Index of Neurons','fontweight','bold');ylabel('Orientation Selectivity Index of Blood Vessels','fontweight','bold');title('Orientation Selectivity Index','fontweight','bold');ylim([0 1]);xlim([0 1]);

figure();boxplot(sel_selectivity,'Labels',{'Neurons','Vessels'});
set(gca,'fontweight','bold','LineWidth',2)
ylim([0 1]);title('Distribution of Orientation Selectivity Index','fontweight','bold');hold on;
for pi=1:length(idx_pi)
    scatter(1,sel_selectivity(pi,1),'g','LineWidth',2);
    scatter(2,sel_selectivity(pi,2),'r','LineWidth',2);
end
hold off;

% for pi=1:1%length(idx_pi)
%     y=squeeze(lissom.layers{2}.xx8(idx_pi(pi),idx_pj(pi),:));
%     y1=transpose(y);
%     y1(1,9)=y(1,1);
%     t=[45 90 135 180 225 270 315 360];
%     th=[0,0.785398163397448,1.570796326794897,2.356194490192345,3.141592653589793,3.926990816987241,4.712388980384690,5.497787143782138,6.283185307179586];
%     figure();subplot(211);mmpolar(th,y1,'g-o'),title('Neural Layer ');
%     set(gca,'fontweight','bold','LineWidth',2)
%     y=squeeze(lissom.layers{3}.xx8(idx_pi(pi),idx_pj(pi),:));
%     y1=transpose(y);
%     y1(1,9)=y(1,1);
%     subplot(212);mmpolar(th,y1,'r-o'),title('Vascular layer');%(strcat('layer 2',{' '},num2str(idx_pi(pi)),num2str(idx_pj(pi))));
%     set(gca,'fontweight','bold','LineWidth',2)
%     %     filename = ['filename' num2str(idx_pi(pi)),'  ',  num2str(idx_pj(pi)) '.jpg'];
%     %     saveas(fig(idx_pi(pi),idx_pj(pi)),filename)    % here you save the figure
% end

figure();
for pi=1:length(idx_pi)
    scatter(lissom.layers{2}.diavg(idx_pi(pi),idx_pj(pi)),lissom.layers{3}.diavg(idx_pi(pi),idx_pj(pi)),'b','LineWidth',2);hold on;
    sel_di(pi,1)=lissom.layers{2}.diavg(idx_pi(pi),idx_pj(pi));
    sel_di(pi,2)=lissom.layers{3}.diavg(idx_pi(pi),idx_pj(pi));
end
plot(0:0.1:1,0:.1:1,'r','LineWidth',2);
set(gca,'fontweight','bold','LineWidth',2)
xlabel('Directionality Index of Neurons','fontweight','bold');ylabel('Directionality Index of Blood Vessels','fontweight','bold');title('Directionality Index','fontweight','bold');ylim([0 1]);xlim([0 1]);

figure();boxplot(sel_di,'Labels',{'Neurons','Vessels'});
set(gca,'fontweight','bold','LineWidth',2)
ylim([0 1]);title('Distribution of Directionality Index');hold on;
for pi=1:length(idx_pi)
    scatter(1,sel_di(pi,1),'g','LineWidth',2);
    scatter(2,sel_di(pi,2),'r','LineWidth',2);
end
hold off;

figure();
for i=1:size(X1,2)
    plot(lissom.layers{2, 1}.converg(i,2:end));hold on;
end
end