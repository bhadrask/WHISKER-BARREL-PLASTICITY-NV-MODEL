function [yy]=inst_response(lissom,opts,X1,Z_constraint,ATP,Amax)
yy=[];
ATP1=reshape(ATP,lissom.layers{2}.dim(1), lissom.layers{2}.dim(2));
inpsel=X1;

for jj=2:numel(lissom.layers)
    [lissom] = rf_m(lissom,double(inpsel),jj,opts);
    [lissom] = activate(lissom,jj,opts);
    
    [lissom] = lat_dynamics_final_vas_intr(lissom,jj,opts,Z_constraint,ATP1,Amax);
    
    yy=lissom.layers{jj}.Zold;
    
    
end