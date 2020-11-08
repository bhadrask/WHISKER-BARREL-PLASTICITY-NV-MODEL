function [proceed] = check(lissom,a)
chx=[];chy=[];           %stride=3;
proceed=[1 1];


SL1=lissom.layers{2}.dim(1);
SL2=lissom.layers{3}.dim(1);
RFL1=lissom.layers{2}.rf(1);
RFL2=lissom.layers{3}.rf(1);
PF=lissom.layers{2}.rf_fb(1);
val1=min(RFL2,PF);
cond1=ceil((size(a{1},1)-RFL1)/lissom.layers{2}.stride+1);
cond2=ceil((SL1-RFL2)/lissom.layers{3}.stride+1);
cond3=ceil((SL1-PF)/lissom.layers{3}.stride+1);
 if SL1<cond1 || SL2<max(cond2,cond3)
     proceed=0;
     error('Size do not match');
 end


% for ii = 2 : numel(lissom.layers)% for each layer
%     if ii==3
%         chx(ii-1)=lissom.layers{ii}.dim(1)-(ceil((lissom.layers{ii-1}.dim(1)-lissom.layers{ii}.rf(1))/lissom.layers{ii}.stride)+1);  %10-((ceil(32-6))/3)+1 == 0
%         chy(ii-1)=lissom.layers{ii}.dim(2)-(ceil((lissom.layers{ii-1}.dim(2)-lissom.layers{ii}.rf(2))/lissom.layers{ii}.stride)+1);
%     end
%     if ii==2
%         chx(ii-1)=lissom.layers{ii}.dim(1)-(ceil((size(a{1},1)-lissom.layers{ii}.rf(1))/lissom.layers{ii}.stride)+1);  %10-((ceil(32-6))/3)+1 == 0
%         chy(ii-1)=lissom.layers{ii}.dim(2)-(ceil((size(a{1},2)-lissom.layers{ii}.rf(2))/lissom.layers{ii}.stride)+1);  %10-((ceil(32-6))/3)+1 == 0
%     end
% end
% for j=1:numel(lissom.layers)-1 % n2=(n1-rf)/stride + 1
%     if chx(j)<0||chy(j)<0    % 0||0 = 0 so the proceed will not change to 0 // proceed = 1
%         proceed(1,j)=0; %if both are not 0, condition wil not satisfy and proceed = 0
%     end
% end
end