function [resp]=gactive(state,alphal,alphau)
ind1=find(state<=alphal);
ind2=find(state>=alphau);
if numel(ind1)~=0
state(ind1)=0;
end
if numel(ind2)~=0
state(ind2)=1;
end
ind3=find(state~=0 & state~=1);
if numel(ind3)~=0
state(ind3)=(state(ind3)-alphal)/(alphau-alphal);
end
resp=state;
end