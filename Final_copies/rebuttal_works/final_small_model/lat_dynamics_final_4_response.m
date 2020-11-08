function [lissom] =  lat_dynamics_final_4_response(lissom,ii,opts,input,ATP,ATPmax.thr)
lissom_aff = lissom.layers{ii}.waff*input;ttt=0;
for k= 1:1  %opts.niter = 8 %number of iterations for zo to adapt its response for 1 image input
    positive=lissom.layers{ii}.wexc*lissom.layers{ii}.Zold;
    negative=lissom.layers{ii}.winh*lissom.layers{ii}.Zold;
 
%     lissom.layers{ii}.Z  =opts.p(1,ii-1)*lissom_aff +(opts.q(1,ii-1)*positive)-(opts.r(1,ii-1)*negative);  %response total=aff + exci + inhi
% thr=2.1;%2.46;

threshold= fnctn(thr,ATP);
lissom.layers{ii}.Zold  = opts.p*lissom_aff +(opts.q*positive)-(opts.r*negative)-threshold;
  
    hi=ATPmax;%opts.p+opts.q-opts.r;
 lo=hi/6;
%       
 lissom.layers{ii}.Zold_passed  = gactive(lissom.layers{ii}.Zold ,lo,hi);
disp(['maximum value before threshold = ',num2str(max(max(lissom.layers{ii}.Zold)))]);

disp(['maximum neural value = ',num2str(max(max(lissom.layers{ii}.Zold_passed)))]);

end
%      figure(11);subplot(1,2,2);imagesc(reshape(lissom.layers{ii}.Zold,lissom.layers{ii}.dim));title(['LEVEL ',num2str(ii)]);colorbar;pause(0.01);
end