function [lissom] = activate(lissom,ii,opts)
for i=1: lissom.layers{ii}.dim(1)
    for j=1: lissom.layers{ii}.dim(2)
%           if numel(find(abs(lissom.layers{ii}.X1{i,j})>0))~=0
%         lissom.layers{ii}.X1{i,j}=normalize(lissom.layers{ii}.X1{i,j});
%           end
   lissom.layers{ii}.Z(i,j) = opts.p*sum(sum(lissom.layers{ii}.X1{i,j}.* lissom.layers{ii}.waff{i,j}));      
%     norm(lissom.layers{ii}.waff{i,j})
%     norm(lissom.layers{ii}.X1{i,j})
    end
end
end
%lissom.layers{ii}.Z = opts.p*sum(sum(repmat(lissom.layers{ii}.X,10,10).*cell2mat(lissom.layers{ii}.waff)));
 