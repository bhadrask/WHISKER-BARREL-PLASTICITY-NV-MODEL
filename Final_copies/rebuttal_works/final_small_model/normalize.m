function [w] =normalize(w,p)
if p==2
    w=w/norm(w);
elseif p==1
    if sum(sum(abs(w)))==0
        w=w/1;
    else
        w=w/sum(sum(abs(w)));
    end
else
    disp('enter valid norm 1 or 2');
end
end