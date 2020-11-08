function [lissom] = feedback_new(lissom,ii,opts)
rf1=lissom.layers{ii}.rf_fb(1)-1;
rf_left=floor(rf1/2);
rf_right=ceil(rf1/2);
cx=lissom.layers{2}.vessel_adjust_indx.cx;
cy=lissom.layers{2}.vessel_adjust_indx.cy;
G1=gaussian_filter_wt(lissom.layers{ii}.rf_fb(1),ceil(lissom.layers{ii}.rf_fb(1)/2));

Mf=zeros(lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
Basemat=ones(lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
Z2d=reshape(lissom.layers{ii+1}.Zold,lissom.layers{ii+1}.dim(1),lissom.layers{ii+1}.dim(2));
cnt=0;
for j=1: lissom.layers{ii+1}.dim(1)
    for i=1: lissom.layers{ii+1}.dim(2)
        cnt=cnt+1;
        G2=zeros(lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
        whole_int_y=cy(cnt)-rf_left:cy(cnt)+rf_right;
        whole_int_x=cx(cnt)-rf_left:cx(cnt)+rf_right;
        final_y=whole_int_y(whole_int_y>0);
        final_x=whole_int_x(whole_int_x>0);
        ynum=abs(numel(final_y)-size(G1,1))+1;
        xnum=abs(numel(final_x)-size(G1,2))+1;
       G2(final_y,final_x)=G1(ynum:end,xnum:end);
         G=G2(1:lissom.layers{ii}.dim(1),1:lissom.layers{ii}.dim(1)); 
        v=Z2d(cnt);
%      figure(1); imagesc(G2);pause(0.1);
        X=v*Basemat;%repmat(v,lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
        Mf = Mf + opts.s*X.*G;
%         figure(1);imagesc(G); pause(0.3);
       
    end
end
lissom.layers{ii}.Zfb=reshape(Mf,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
%  figure(1);imagesc(Mf); pause(0.3);
end