function [lissom] = make_unit_norm(lissom,ii,p)
 for j = 1 : lissom.layers{ii}.dim(1) %each row in layer
         for k = 1 : lissom.layers{ii}.dim(2) %each column in layer
lissom.layers{ii}.waff{j,k}=normalize((lissom.layers{ii}.waff{j,k}),p);

lissom.layers{ii}.wnbrexc{j,k}=normalize((lissom.layers{ii}.wnbrexc{j,k}),p);
lissom.layers{ii}.wnbrinhb{j,k}=normalize((lissom.layers{ii}.wnbrinhb{j,k}),p);
         end
 end