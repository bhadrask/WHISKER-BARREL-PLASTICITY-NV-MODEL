function [lissom] = initialise4usage(lissom,filename)
load(filename);
lissom.layers{2}.waff=lissomaff;
lissom.layers{2}.wnbrexc=lissomexc;
lissom.layers{2}.wnbrinhb=lissominhb;
opts.p=gama1;opts.q=gama2;opts.r=gama3;

 for ii = 2 : numel(lissom.layers)% for each layer
     exc_mat=zeros(lissom.layers{ii}.exe_rad,lissom.layers{ii}.exe_rad);
     inh_mat=zeros(lissom.layers{ii}.inhb_rad,lissom.layers{ii}.inhb_rad);
   if strcmp(lissom.layers{ii}.type, 'v')
       
        % initialize weights for each neuron wrt corresponding rf as waff
        % construct lateral neighbouring matrices and corresponding weight matrices
        
        [xgrid, ygrid] = meshgrid(1:lissom.layers{ii}.dim(2), 1:lissom.layers{ii}.dim(1));  
        
         
       for j = 1 : lissom.layers{ii}.dim(1) %each row in layer
         for k = 1 : lissom.layers{ii}.dim(2) %each column in layer
            
           
%         
          % Make a logical circular region set to 1, the rest to zero
            x = xgrid - k;    % offset the current position
            y = ygrid - j;
            lissom.layers{ii}.nbrexc{j,k} = x.^2 + y.^2 <= lissom.layers{ii}.exe_rad.^2;
            lissom.layers{ii}.outer_inhb{j,k} = x.^2 + y.^2 <= lissom.layers{ii}.inhb_rad.^2;
          %inhibitory response with in the specified radius
           lissom.layers{ii}.nbrinhb{j,k}= ~lissom.layers{ii}.nbrexc{j,k}.*lissom.layers{ii}.outer_inhb{j,k};
           lissom.layers{ii}.inhb_idx{j,k}=find((lissom.layers{ii}.nbrinhb{j,k})>0);
           lissom.layers{ii}.exc_idx{j,k}=find((lissom.layers{ii}.nbrexc{j,k})>0);
          
         end 
       end 
            lissom.layers{ii}.PreSynLat = zeros(lissom.layers{ii}.dim);
       lissom.layers{ii}.Zold =zeros(lissom.layers{ii}.dim);
       lissom.layers{ii}.Z  =zeros(lissom.layers{ii}.dim);
   end
end