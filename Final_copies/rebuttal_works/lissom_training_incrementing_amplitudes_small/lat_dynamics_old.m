function [lissom] = lat_dynamics_old(lissom,ii,opts,Z_constraint)
%number of time steps for setteling is predefined. Instead can be used threshold  
lissom_aff = (lissom.layers{ii}.Z).*(Z_constraint);ttt=0;
for k= 1:opts.niter
   for m=1: lissom.layers{ii}.dim(1)
     for n=1: lissom.layers{ii}.dim(2)
         
     
        positive(m,n) = sum(sum(lissom.layers{ii}.Zold.*lissom.layers{ii}.wnbrexc{m,n}));
        negative(m,n) = sum(sum(lissom.layers{ii}.Zold.*lissom.layers{ii}.wnbrinhb{m,n}));
    end; 
   end;
        lissom.layers{ii}.Z  = lissom_aff +(opts.q*positive)-(opts.r*negative);
  
        if ii>1 && ii<numel(lissom.layers) && opts.control==1
        lissom_fb=lissom.layers{ii}.Zfb;
        lissom.layers{ii}.Z=lissom.layers{ii}.Z+lissom_fb;
        end
%         
  hi=max(max(lissom.layers{ii}.Z));       
%  lo=min(min(lissom.layers{ii}.Z));
%  temp=(hi-lo)*0.5;
%   lo=lo+temp;
%         
 %hi=opts.p+opts.q-opts.r;
 lo=0.6;
%       
 lissom.layers{ii}.Z1  = gactive(lissom.layers{ii}.Z ,lo,hi);
  lissom.layers{ii}.Zold = lissom.layers{ii}.Z1 .*Z_constraint;
 if numel(find(abs(lissom.layers{ii}.Zold)>0))~=0
         lissom.layers{ii}.Zold=normalize(lissom.layers{ii}.Zold);
 end
 delz(k)=norm(lissom.layers{ii}.Zold -ttt);
 ttt=lissom.layers{ii}.Zold;
 if k==opts.niter-1
   lissom.layers{ii}.Zprev=lissom.layers{ii}.Zold;
 end     
        
     
end
%  figure(1)
%  subplot(3,3,ii);title(strcat('layer',num2str(ii)))
%  imagesc(lissom.layers{ii}.Zold); pause(0.001);
  
