function [lissom] = RF_master_new(lissom,ii,aff_sheet1d)
    P=lissom.layers{ii}.zeropad;
    aff_sheet2d=reshape(aff_sheet1d,sqrt(numel(aff_sheet1d)),sqrt(numel(aff_sheet1d)));
    aff_sheet=padarray(aff_sheet2d,[P(1) P(2)]);
    aff_padded=reshape(aff_sheet,size(aff_sheet,1)*size(aff_sheet,2),1);
    RF_master=bsxfun(@times, lissom.layers{ii}.I_aff_master, aff_padded');
    lissom.layers{ii}.RF_master= sparse(RF_master);
    lissom.layers{ii}.padded_input=aff_padded;
end