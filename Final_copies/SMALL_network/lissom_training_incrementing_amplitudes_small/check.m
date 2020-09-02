function [proceed] = check(lissom,a)
chx=[];
chy=[];
%stride=3;
stride1=lissom.layers{2}.stride;
proceed=1;
chx(1)=lissom.layers{2}.dim(1)-(ceil((size(a,1)-lissom.layers{2}.rf(1))/stride1)+1);
chy(1)=lissom.layers{2}.dim(2)-(ceil((size(a,2)-lissom.layers{2}.rf(2))/stride1)+1);
for i=3:numel(lissom.layers)
    stride2=lissom.layers{i}.stride;
    chx(i-1)=lissom.layers{i}.dim(1)-(ceil((lissom.layers{i-1}.dim(1)-lissom.layers{i}.rf(1))/stride2)+1);
    chy(i-1)=lissom.layers{i}.dim(1)-(ceil((lissom.layers{i-1}.dim(2)-lissom.layers{i}.rf(2))/stride2)+1);
end
for j=1:numel(lissom.layers)-1
    if chx(j)<0||chy(j)<0
        proceed=0;
    end
end
