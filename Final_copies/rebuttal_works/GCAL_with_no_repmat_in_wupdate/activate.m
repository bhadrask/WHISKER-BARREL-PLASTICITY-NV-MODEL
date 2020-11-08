function [lissom] = activate(lissom,ii,input)
% for i=1: lissom.layers{ii}.dim(1)       %1 - 10
%     for j=1: lissom.layers{ii}.dim(2)   %1 - 10
%         lissom.layers{ii}.Z(i,j) = opts.p(1,ii-1)*sum(sum(lissom.layers{ii}.X1{i,j}.* lissom.layers{ii}.waff{i,j}));      
%     end  %response of neuron in V1 is product of p(weight of afferent connection) and 
% end  %sum of the product of input and waff(weight)
 lissom.layers{ii}.Z=lissom.layers{ii}.waff*input;
end
 