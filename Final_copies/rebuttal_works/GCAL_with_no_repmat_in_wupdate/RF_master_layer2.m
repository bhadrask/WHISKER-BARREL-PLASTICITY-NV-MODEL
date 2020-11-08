function [lissom] = RF_master_layer2(lissom,ii,opts)
        aff_sheet=lissom.layers{ii-1}.Zold;
        aff_sheet=padarray(aff_sheet,[pad1 pad2]); %[1 1]
       [lissom]=prefix_mats(lissom,aff_sheet,ii,opts);
        RF_master=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),numel(aff_sheet));
   for ii1=1:size(lissom.layers{ii}.I,1)
        RF_master(ii1,lissom.layers{ii}.I(ii1,:))=aff_sheet(lissom.layers{ii}.I(ii1,:));
   end
lissom.layers{ii}.RF_master(ii)= RF_master;
