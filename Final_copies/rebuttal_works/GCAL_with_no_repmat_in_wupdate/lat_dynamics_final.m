function [lissom] =  lat_dynamics_final(lissom,ii,opts,input,ATP,ATPmax)
lissom_aff = lissom.layers{ii}.waff*input;ttt=0;
for k= 1:opts.niter  %opts.niter = 8 %number of iterations for zo to adapt its response for 1 image input
    positive=lissom.layers{ii}.wexc*lissom.layers{ii}.Zold;
    negative=lissom.layers{ii}.winh*lissom.layers{ii}.Zold;
 
%     lissom.layers{ii}.Z  =opts.p(1,ii-1)*lissom_aff +(opts.q(1,ii-1)*positive)-(opts.r(1,ii-1)*negative);  %response total=aff + exci + inhi
thr=2.1;%2.46;

threshold= fnctn(thr,ATP);
lissom.layers{ii}.Z  = opts.p*lissom_aff +(opts.q*positive)-(opts.r*negative)-threshold;
  
    hi=ATPmax;%opts.p+opts.q-opts.r;
 lo=hi/6;
%       
 lissom.layers{ii}.Zold  = gactive(lissom.layers{ii}.Z ,lo,hi);

 delz(k)=norm(lissom.layers{ii}.Zold -ttt);
 ttt=lissom.layers{ii}.Zold;
 if k==opts.niter-1
   lissom.layers{ii}.Zprev=lissom.layers{ii}.Zold;
 end     
%       figure(11);imagesc(reshape(lissom.layers{ii}.Z,lissom.layers{ii}.dim));title(['LEVEL ',num2str(ii)]);colorbar;pause(0.01);

%      figure(11);imagesc(reshape(lissom.layers{ii}.Zold,lissom.layers{ii}.dim));title(['LEVEL ',num2str(ii)]);colorbar;pause(0.03);
end
%      figure(11);subplot(1,2,2);imagesc(reshape(lissom.layers{ii}.Zold,lissom.layers{ii}.dim));title(['LEVEL ',num2str(ii)]);colorbar;pause(0.01);
end