function [lissom] = RF_master(lissom,ii,aff_sheet1d)
P=lissom.layers{ii}.zeropad;
aff_sheet2d=reshape(aff_sheet1d,sqrt(numel(aff_sheet1d)),sqrt(numel(aff_sheet1d)));
        aff_sheet=padarray(aff_sheet2d,[P(1) P(2)]); 
        aff_padded=reshape(aff_sheet,size(aff_sheet,1)*size(aff_sheet,2),1);
        Temp_pad=repmat(aff_padded',lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),1);
        RF_master=zeros(lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),numel(aff_sheet));
   for ii1=1:size(lissom.layers{ii}.I,1)
        RF_master(ii1,lissom.layers{ii}.I(ii1,:))=aff_sheet(lissom.layers{ii}.I(ii1,:));
   end
lissom.layers{ii}.RF_master= sparse(RF_master);

lissom.layers{ii}.padded_input=aff_padded;