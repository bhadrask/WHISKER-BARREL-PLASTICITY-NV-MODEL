function [lissom] = initialise(lissom,xdim,c,opts)

 for ii = 2 : numel(lissom.layers)% for each layer
   if strcmp(lissom.layers{ii}.type, 'v')
       
        % initialize weights for each neuron wrt corresponding rf as waff
        % construct lateral neighbouring matrices and corresponding weight matrices
        
        [xgrid, ygrid] = meshgrid(1:lissom.layers{ii}.dim(2), 1:lissom.layers{ii}.dim(1));  
        
         
       for j = 1 : lissom.layers{ii}.dim(1) %each row in layer
         for k = 1 : lissom.layers{ii}.dim(2) %each column in layer
             if ii>1 && ii<numel(lissom.layers) && opts.control==1
              lissom.layers{ii}.wfdbk{j,k}=1e-3*normalize(rand(lissom.layers{ii+1}.rf));
             end
           lissom.layers{ii}.waff{j,k} =(rand(lissom.layers{ii}.rf));
           
           wexc{j,k} = (rand(lissom.layers{ii}.dim)); 
           winhib{j,k} = (rand(lissom.layers{ii}.dim)); 
           
          % Make a logical circular region set to 1, the rest to zero
            x = xgrid - k;    % offset the current position
            y = ygrid - j;
            lissom.layers{ii}.nbrexc{j,k} = x.^2 + y.^2 <= lissom.layers{ii}.exe_rad.^2;
          %inhibitory response with in the specified radius
           lissom.layers{ii}.nbrinhb{j,k}= ~lissom.layers{ii}.nbrexc{j,k};%global inhibitory response
           lissom.layers{ii}.wnbrexc{j,k} = lissom.layers{ii}.nbrexc{j,k}.*  wexc{j,k};
           lissom.layers{ii}.wnbrinhb{j,k} = lissom.layers{ii}.nbrinhb{j,k}.* winhib{j,k};
         end 
       end 
       lissom.layers{ii}.PreSynLat = zeros(lissom.layers{ii}.dim);
       lissom.layers{ii}.Zold =zeros(lissom.layers{ii}.dim);
       lissom.layers{ii}.Z  =zeros(lissom.layers{ii}.dim);
       lissom.layers{ii}.Zabove=zeros(lissom.layers{ii}.dim);
   end
end
     
      