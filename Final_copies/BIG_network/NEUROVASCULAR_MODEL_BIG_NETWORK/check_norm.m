function [x,y,z,q] = check_norm(lissom,ii,p,opts)
x=0;y=0;z=0;q=0;
for j = 1 : lissom.layers{ii}.dim(1) %each row in layer
        for k = 1 : lissom.layers{ii}.dim(2) %each column in layer
waff1=sum(sum((lissom.layers{ii}.waff{j,k}),p));

wexc1=sum(sum((lissom.layers{ii}.wnbrexc{j,k}),p));
winh1=sum(sum((lissom.layers{ii}.wnbrinhb{j,k}),p));
if ii>1 && ii<numel(lissom.layers) && opts.control==1
             wfdbk1=sum(sum((lissom.layers{ii}.wfdbk{j,k}),p));
 if (wfdbk1-1)>1e-5
   q=q+1;
end  
end
if (waff1-1)>1e-5
   x=x+1;
end
if (wexc1-1)>1e-5
   y=y+1;
end
if (winh1-1)>1e-5
   z=z+1;
end
        end
end