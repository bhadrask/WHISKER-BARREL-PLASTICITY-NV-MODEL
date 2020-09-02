function [yy]=inst_response(lissom,opts,X1,Z_constraint,label,ATP,ATPmax)
yy=[];xx=[];
label(length(label)+1)=-1;
for i=1:size(X1,3)
    inpsel=X1(:,:,i);
    figure(2); subplot(2,1,1);imagesc(inpsel);
    for jj=2:numel(lissom.layers)
        [lissom] = rf_m(lissom,double(inpsel),jj,opts);
        [lissom] = activate(lissom,jj,opts);
        
        [lissom]=lat_dynamics_final_4_response(lissom,jj,opts,Z_constraint,ATP,ATPmax);
        
        yy=lissom.layers{jj}.Zold;
        
        xx=cat(3,xx,lissom.layers{jj}.Zold);
        
    end
    figure(2);
    subplot(2,1,2);imagesc(yy);colorbar;caxis([0 1]);title(['input',num2str(i),' max=',num2str(max(max(yy)))]); pause(0.0001)
    for kk=2:numel(lissom.layers)
            lissom.layers{kk}.Zold = zeros(lissom.layers{kk}.dim);
            lissom.layers{kk}.Z = zeros(lissom.layers{kk}.dim);% reset activation to zero before new input
    end
    
end;
[gg,q]=max(xx,[],3);
% figure; imagesc(labeldir(q));title('layer2 direction label');colorbar;
loc= gg<=0;
q(loc)=length(label);
figure(); imagesc(label(q));title('barrels');colorbar;