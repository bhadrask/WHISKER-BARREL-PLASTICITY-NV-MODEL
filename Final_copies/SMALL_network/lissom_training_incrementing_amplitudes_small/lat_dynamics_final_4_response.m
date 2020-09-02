function [lissom] = lat_dynamics_final_4_response(lissom,ii,opts,ATP,ATPmax)
%number of time steps for setteling is predefined. Instead can be used threshold  
lissom_aff = (lissom.layers{ii}.Z);ttt=0;
for k= 1:1
%     lissom.layers{ii}.Z2=repmat(lissom.layers{ii}.Zold,1,1,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2));

   for m=1: lissom.layers{ii}.dim(1)
     for n=1: lissom.layers{ii}.dim(2)
         zexc=lissom.layers{ii}.Zold(lissom.layers{ii}.exc_idx{m,n});
         zinhb=lissom.layers{ii}.Zold(lissom.layers{ii}.inhb_idx{m,n});
     
        positive(m,n) = sum(sum(zexc.*lissom.layers{ii}.wnbrexc{m,n}));
        negative(m,n) = sum(sum(zinhb.*lissom.layers{ii}.wnbrinhb{m,n}));
    end; 
   end;
   




        lissom.layers{ii}.Z  = lissom_aff +opts.q*(positive)-opts.r*(negative)+ATP;
  
        if ii>1 && ii<numel(lissom.layers) && opts.control==1
        lissom_fb=lissom.layers{ii}.Zfb;
        lissom.layers{ii}.Z=lissom.layers{ii}.Z+lissom_fb;
        end

%         
 hi=ATPmax;
 lo=hi/6;
%    
disp(['maximum value before threshold = ',num2str(max(max(lissom.layers{ii}.Z)))]);
 lissom.layers{ii}.Zold  = gactive(lissom.layers{ii}.Z ,lo,hi);
disp(['maximum neural value = ',num2str(max(max(lissom.layers{ii}.Zold)))]);

   
end

