function [w] =normalize(w,p)
% c = [0.9501    0.4860    0.4565; 0.2311    0.8913    0.0185; 0.6068    0.7621    0.8214]
% w=c;
% s1=size(w);
% wreshape=reshape(w,s1(1)*s1(2),1);
% wsqr=wreshape.^2;
% normsum=sum(wsqr); 
% wnorm=wreshape./sqrt(normsum);
% w=reshape(wnorm,s1(1),s1(2));
if p==2
 w=w/norm(w);
elseif p==1
w=w/sum(sum(abs(w)));
else
    disp('enter valid norm');
end
end