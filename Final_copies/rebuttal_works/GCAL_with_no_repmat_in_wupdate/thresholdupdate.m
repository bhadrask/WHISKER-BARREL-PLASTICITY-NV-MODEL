function [lissom] =thresholdupdate(lissom,ii)
lissom.layers{ii}.oldthreshold = lissom.layers{ii}.threshold;
if ii==2    
    lissom.layers{ii}.ave_res = 0.999*lissom.layers{ii}.prev_ave_res + 0.001*lissom.layers{ii}.rectify_res;
    lissom.layers{ii}.threshold= lissom.layers{ii}.oldthreshold+ 0.0001*lissom.layers{ii}.ave_res - (0.0001*0.024);
    lissom.layers{ii}.prev_ave_res= lissom.layers{ii}.ave_res;    
end
if ii==3    
    lissom.layers{ii}.ave_res= 0.999*lissom.layers{ii}.prev_ave_res+ 0.001*lissom.layers{ii}.rectify_res;
    lissom.layers{ii}.threshold= lissom.layers{ii}.oldthreshold+ 0.0001*lissom.layers{ii}.ave_res - (0.0001*0.024);
    lissom.layers{ii}.prev_ave_res= lissom.layers{ii}.ave_res;    
end
end