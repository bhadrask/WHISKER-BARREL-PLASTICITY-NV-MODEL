function [lissom] = lat_dynamics_final(lissom,ii,opts,ATP,ATPmax)
%number of time steps for setteling is predefined. Instead can be used threshold  
lissom_aff = (lissom.layers{ii}.Z);ttt=0;
for k= 1:opts.niter
%     lissom.layers{ii}.Z2=repmat(lissom.layers{ii}.Zold,1,1,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));

   for m=1: lissom.layers{ii}.dim(1)
     for n=1: lissom.layers{ii}.dim(2)
         zexc=lissom.layers{ii}.Zold(lissom.layers{ii}.exc_idx{m,n});
         zinhb=lissom.layers{ii}.Zold(lissom.layers{ii}.inhb_idx{m,n});
     
        positive(m,n) = sum(sum(zexc.*lissom.layers{ii}.wnbrexc{m,n}));
        negative(m,n) = sum(sum(zinhb.*lissom.layers{ii}.wnbrinhb{m,n}));
    end; 
   end;
   

%    Zexc= cellfun(@(x) lissom.layers{ii}.Zold(x),lissom.layers{ii}.exc_idx, 'UniformOutput',false);
%      Zinhb= cellfun(@(x) lissom.layers{ii}.Zold(x),lissom.layers{ii}.inhb_idx, 'UniformOutput',false);
%      positive=cell2mat(cellfun(@(x,y) sum(sum(x.*y)), Zexc,lissom.layers{ii}.wnbrexc, 'UniformOutput',false));
%      negative=cell2mat(cellfun(@(x,y) sum(sum(x.*y)), Zinhb,lissom.layers{ii}.wnbrinhb, 'UniformOutput',false));
  thr=2.1;%2.46;

threshold= fnctn(thr,ATP);
lissom.layers{ii}.Z  = opts.p*lissom_aff +(opts.q*positive)-(opts.r*negative)-threshold;
%         lissom.layers{ii}.Z  = lissom_aff +opts.q*(positive)-opts.r*(negative)+ATP;
  
        if ii>1 && ii<numel(lissom.layers) && opts.control==1
        lissom_fb=lissom.layers{ii}.Zfb;
        lissom.layers{ii}.Z=lissom.layers{ii}.Z+lissom_fb;
        end
%         
%  hi=max(max(lissom.layers{ii}.Z));       
%  lo=min(min(lissom.layers{ii}.Z));
%  temp=(hi-lo)*0.5;
%   lo=lo+temp;
%         
 hi=ATPmax;%opts.p+opts.q-opts.r;
 lo=hi/6;
%       
 lissom.layers{ii}.Zold  = gactive(lissom.layers{ii}.Z ,lo,hi);
%  if numel(find(abs(lissom.layers{ii}.Zold)>0))~=0
%          lissom.layers{ii}.Zold=normalize(lissom.layers{ii}.Zold);
%  end
 delz(k)=norm(lissom.layers{ii}.Zold -ttt);
 ttt=lissom.layers{ii}.Zold;
 if k==opts.niter-1
   lissom.layers{ii}.Zprev=lissom.layers{ii}.Zold;
 end     
%           figure(2)
%    imagesc(lissom.layers{ii}.Zold) ;title(num2str(k)) ;pause(0.1);

   
end
%    figure(3);
%    plot(delz)
%  figure(1)
%  subplot(3,3,ii);title(strcat('layer',num2str(ii)))
%  imagesc(lissom.layers{ii}.Zold); %pause(0.1);
  
