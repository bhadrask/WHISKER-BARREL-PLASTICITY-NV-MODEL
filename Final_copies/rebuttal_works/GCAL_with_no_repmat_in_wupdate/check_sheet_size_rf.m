function [proceed] = check_sheet_size_rf(lissom,a)
chx=[];chy=[];           %stride=3;
proceed=[1 1];


SL1=lissom.layers{2}.dim(1);

RFL1=lissom.layers{2}.rf(1);
cond1=ceil((size(a,1)-RFL1)/lissom.layers{2}.stride+1);
if SL1<cond1
    proceed(1)=0;
    error(['SL1 should atleast be ',num2str(cond1)]);
    elseif (SL1-cond1)>5
    error(['too much zero padding in input layer, reduce SL1 to ',num2str(cond1)]);
end

if numel(lissom.layers)>2
SL2=lissom.layers{3}.dim(1);
RFL2=lissom.layers{3}.rf(1);
cond2=ceil((SL1-RFL2)/lissom.layers{3}.stride+1);

if SL2<cond2
      proceed(2)=0;
    error(['SL2 should atleast be ',num2str(cond2)])

    elseif (SL2-cond2)>3
    error(['too much zero padding in neural layer, reduce SL2 to ',num2str(cond2)]);
end

end









end