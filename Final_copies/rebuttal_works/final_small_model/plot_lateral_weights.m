function [lissom]=plot_lateral_weights(lissom)
for ii=2:numel(lissom.layers)
    locstot=[];
    locsmegaexch=[];
    locsmegainh=[];
    for m=1: lissom.layers{ii}.dim(1)      %1 to 10
        for n=1: lissom.layers{ii}.dim(2)
            locs=zeros(lissom.layers{ii}.dim);
            locs1=zeros(lissom.layers{ii}.dim);
            zexc=(lissom.layers{ii}.exc_idx{m,n});    %z values in the excitatory radius
            zinhb=(lissom.layers{ii}.inhb_idx{m,n});
            locs(zexc)=lissom.layers{ii}.wnbrexc{m,n};
            locs1(zinhb)=lissom.layers{ii}.wnbrinhb{m,n};
            total=locs-locs1;
            locstot{m,n}=total;
            locsmegaexch{m,n}=locs;
            locsmegainh{m,n}=locs1;
            %              locsmega=cat(2,locsmega,locs);
            %             figure(1);subplot(211);imagesc(locs);title(['layer', num2str(ii),'lateral exc weights']);
            %             figure(1);subplot(212);imagesc(locs1);title(['layer', num2str(ii),'lateral inh weights']);pause(0.1)
        end
    end
    lissom.layers{ii}.storewlateral_exch=cat(3,lissom.layers{ii}.storewlateral_exch,locsmegaexch);
    lissom.layers{ii}.storewlateral_inhi=cat(3,lissom.layers{ii}.storewlateral_inhi,locsmegainh);
    lissom.layers{ii}.storewaff=cat(3,lissom.layers{ii}.storewaff,lissom.layers{ii, 1}.waff);    
     lissom.layers{ii}.storethr=cat(3,lissom.layers{ii}.storethr,lissom.layers{ii, 1}.threshold);    
%    figure(ii);subplot(131);imagesc(cell2mat(locsmegaexch));title(['layer', num2str(ii),'lateral exc weights']);
% figure(ii);subplot(132);imagesc(cell2mat(locsmegainh));title(['layer', num2str(ii),'lateral inh weights']);
% figure(ii);subplot(133);imagesc(cell2mat(lissom.layers{ii, 1}.waff)); title(['layer', num2str(ii),'aff']);
end
