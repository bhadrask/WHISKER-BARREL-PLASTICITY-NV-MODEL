function [lissom] = find_neural_avg(lissom,ii,type)
pixel_spacing=3000/lissom.layers{2}.dim(1);%3000 micrometer sheet size
rf1=round(600/pixel_spacing)-1;% considering a sheet size of 3mmx3mm, neural avg is found around 600micromx600microm
rf_left=floor(rf1/2);
rf_right=ceil(rf1/2);
cx=lissom.layers{2}.vessel_adjust_indx.cx;
cy=lissom.layers{2}.vessel_adjust_indx.cy;
G1=ones(rf1+1,rf1+1);

if type==1
    Z_osi=lissom.layers{ii}.Selectivity;
    Z_op=lissom.layers{ii}.sigmap;
    Z_DI=lissom.layers{ii}.di;
    Zxxnorm=lissom.layers{ii}.xxnorm;
elseif type==2
    Zxxnorm=lissom.layers{ii}.xxnorm;
end
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
        
        %      figure(1); imagesc(G2);pause(0.1
        if type==1
            X_osi(i,j)=sum(sum(Z_osi.*G))/numel(find(G>0));%repmat(v,lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
            X_op(i,j)=sum(sum(Z_op.*G))/numel(find(G>0));
            X_DI(i,j)=sum(sum(Z_DI.*G))/numel(find(G>0));
            %         figure(1);imagesc(G); pause(0.3);
            G3d=repmat(G,1,1,size(Zxxnorm,3));
            sw=Zxxnorm.*G3d;
            sw1=sum(sum(sw))/numel(find(G>0));
            YY(i,j,:)=squeeze(sw1);
        elseif type==2
            G3d=repmat(G,1,1,size(Zxxnorm,3));
            sw=Zxxnorm.*G3d;
            sw1=sum(sum(sw))/numel(find(G>0));
            YY(i,j,:)=squeeze(sw1);
        end
    end
end
if type==1
    lissom.layers{ii}.Avg_neural_OSI=X_osi;
    lissom.layers{ii}.Avg_neural_DI=X_DI;
    lissom.layers{ii}.Avg_neural_OP=X_op;
    lissom.layers{ii}.Avg_neural_sheet_resp=YY;
elseif type==2
    lissom.layers{ii}.Avg_neural_sheet_resp=YY;
end
%  figure(1);imagesc(Mf); pause(0.3);
end