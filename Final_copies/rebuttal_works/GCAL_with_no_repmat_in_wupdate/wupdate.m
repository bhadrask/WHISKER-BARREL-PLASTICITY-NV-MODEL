function [lissom] =wupdate(lissom,ii,opts)
%     Znodeaff=repmat(lissom.layers{ii}.Zold,1,size(lissom.layers{ii}.waff,2));
%     afftemp = lissom.layers{ii}.waff+ opts.etaaff(1,ii-1) * lissom.layers{ii}.RF_master .*Znodeaff;
afftemp = lissom.layers{ii}.waff+ opts.etaaff(1,ii-1) * lissom.layers{ii}.RF_master .*lissom.layers{ii}.Zold;

    Sumnorm=sum(abs(afftemp),2);
    % normmat=repmat(Sumnorm,1,size(lissom.layers{ii}.waff,2));
    normmat= Sumnorm(:,ones(1,size(lissom.layers{ii}.waff,2)));
    lissom.layers{ii}.waff= afftemp./normmat;          %Afferent weight updation
    if ii==2
% Zprev_mat=repmat(lissom.layers{ii}.Zprev',lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),1);
%     EXC_Zpost=Zprev_mat.*lissom.layers{ii}.I_exc_master;
     EXC_Zpost=lissom.layers{ii}.Zprev'.*lissom.layers{ii}.I_exc_master;
    elseif ii==3
%     Zold_mat=repmat(lissom.layers{ii}.Zold',lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(1),1);
%     EXC_Zpost=Zold_mat.*lissom.layers{ii}.I_exc_master;
EXC_Zpost=lissom.layers{ii}.Zold'.*lissom.layers{ii}.I_exc_master;
    end
%     Znode=Znodeaff(:,1:size(lissom.layers{ii}.wexc,2));%repmat(lissom.layers{ii}.Zold,1,size(lissom.layers{ii}.wexc,2));


%     %  EXC_Zpost=sparse(EXC_Zpost);
    exctemp  = lissom.layers{ii}.wexc + opts.etaexc* EXC_Zpost.*lissom.layers{ii}.Zold;
    excnorm=sum(abs(exctemp),2);
%     exnorm_mat=excnorm(:,ones(1,size(lissom.layers{ii}.wexc,2)));%repmat(excnorm,1,size(lissom.layers{ii}.wexc,2));
    lissom.layers{ii}.wexc= exctemp./ excnorm;%exnorm_mat;
    if ii==2
%         INH_Zpost=Zprev_mat.*lissom.layers{ii}.I_inh_master;
       INH_Zpost= lissom.layers{ii}.Zprev'.*lissom.layers{ii}.I_inh_master;
    elseif ii==3
%     INH_Zpost=Zold_mat.*lissom.layers{ii}.I_inh_master;
  INH_Zpost=lissom.layers{ii}.Zold'.*lissom.layers{ii}.I_inh_master;
    end

    %   INH_Zpost=sparse(INH_Zpost);
    inhtemp  = lissom.layers{ii}.winh+ opts.etainh* INH_Zpost.*lissom.layers{ii}.Zold;
    inhnorm=sum(abs(inhtemp),2);
%     inhnorm_mat=inhnorm(:,ones(1,size(lissom.layers{ii}.winh,2)));%repmat(inhnorm,1,size(lissom.layers{ii}.winh,2));
    lissom.layers{ii}.winh= inhtemp./inhnorm;%inhnorm_mat;
end