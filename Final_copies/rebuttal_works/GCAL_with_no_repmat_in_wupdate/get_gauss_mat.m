function [lissom] = get_gauss_mat(lissom,ii)
sigma(1)=lissom.layers{ii-1}.rf_fb(1);
sigma(2)=lissom.layers{ii-1}.rf_fb(2);

[xgrid, ygrid] = meshgrid(1:lissom.layers{ii}.dim(2), 1:lissom.layers{ii}.dim(1));  

cnt=0;
G=zeros(lissom.layers{ii+1}.dim(1)*lissom.layers{ii+1}.dim(2));
for i=1: lissom.layers{ii+1}.dim(1)
    for j=1: lissom.layers{ii+1}.dim(2)
        cnt=cnt+1;
        x = (((xgrid-i) + (ygrid-j))/sigma(1)).^2;
        y = (((xgrid-i) - (ygrid-j))/sigma(2)).^2;
        gg = amp*exp(-(x+y)/2);
         gg(gg<1e-6)=0;
      G(count,:)=reshape(gg,1,lissom.layers{ii+1}.dim(1)*lissom.layers{ii+1}.dim(2));
    end
end
lissom.layers{ii}.Gaussian_mat=G;

end