function [lissom] =wupdate(lissom,ii,opts)
lissom.layers{ii}.oldwaff=lissom.layers{ii}.waff;
  for i=1: lissom.layers{ii}.dim(1)
    for j=1: lissom.layers{ii}.dim(2)
         if numel(find(abs(lissom.layers{ii}.X1{i,j})>0))~=0
        lissom.layers{ii}.X1{i,j}=normalize(lissom.layers{ii}.X1{i,j},1);
          end
    lissom.layers{ii}.PreSynLatexc{i,j} = lissom.layers{ii}.Zold(lissom.layers{ii}.exc_idx{i,j});
    lissom.layers{ii}.PreSynLatinhb{i,j} =lissom.layers{ii}.Zold(lissom.layers{ii}.inhb_idx{i,j});
   
    %Afferent weight updation
    
    afftemp = lissom.layers{ii}.waff{i,j}+opts.etaaff* lissom.layers{ii}.X1{i,j} *lissom.layers{ii}.Zold(i,j);
    lissom.layers{ii}.waff{i,j}= afftemp ;%./(sum(sum(afftemp)));
   
    %lateral weight updation
    exctemp  = lissom.layers{ii}.wnbrexc{i,j}+opts.etaexc*lissom.layers{ii}.PreSynLatexc{i,j}* lissom.layers{ii}.Zold(i,j);
    lissom.layers{ii}.wnbrexc{i,j}= exctemp;%./sum(sum(exctemp));
    
    inhibtemp  = lissom.layers{ii}.wnbrinhb{i,j} +opts.etainhib*lissom.layers{ii}.PreSynLatinhb{i,j}* lissom.layers{ii}.Zold(i,j);
    lissom.layers{ii}.wnbrinhb{i,j}= inhibtemp;%./sum(sum(inhibtemp));
    %feedback weight updation
%     if ii>1 && ii<numel(lissom.layers) && opts.control==1
%        fdbktemp = lissom.layers{ii}.wfdbk{i,j}+opts.etaaff*lissom.layers{ii}.X2{i,j} *lissom.layers{ii}.Zold(i,j);
%        lissom.layers{ii}.wfdbk{i,j}= fdbktemp;%./(sum(sum(fdbktemp)));
%     
%     end
    end;
   end; 
%     figure(1); 
%      %subplot(2,1,1); imagesc(inpsel);pause(0.1);title('input');
% subplot(3,3,ii+2); imagesc(cell2mat(lissom.layers{ii}.waff));title(strcat('w_aff layer',num2str(ii)));%pause(0.1);title('w-aff'); % colormap(gray);
% if opts.control==1
% subplot(3,3,6); imagesc(cell2mat(lissom.layers{2}.wfdbk));title(strcat('w_fdbk layer',num2str(ii)));
% end
%subplot(3,1,3);surf(abs(cell2mat(lissom.layers{ii}.waff)-cell2mat(lissom.layers{ii}.oldwaff)));%pause(0.05)
%     subplot(2,1,2); imagesc(cell2mat(lissom.layers{ii}.wnbrexc));pause(0.1); title('w-exc'); % colormap(gray);
% %     subplot(3,2,4);imagesc(cell2mat(lissom.layers{ii}.wnbrinhb));pause(0.1); title('w-inhb');
%     pause(1);

end
  
    
