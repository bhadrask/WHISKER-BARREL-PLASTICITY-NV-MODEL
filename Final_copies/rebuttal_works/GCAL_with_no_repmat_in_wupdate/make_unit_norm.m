function [lissom] = make_unit_norm(lissom,ii)
for j = 1 : lissom.layers{ii}.dim(1) %each row in layer
    for k = 1 : lissom.layers{ii}.dim(2) %each column in layer
            lissom.layers{ii}.waff{j,k}=normalize((lissom.layers{ii}.waff{j,k}),1);
            lissom.layers{ii}.wnbrexc{j,k}=normalize((lissom.layers{ii}.wnbrexc{j,k}),1);
            lissom.layers{ii}.wnbrinhb{j,k}=normalize((lissom.layers{ii}.wnbrinhb{j,k}),1);
    end  %normalize weights
end