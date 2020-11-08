function plotcaps(f,n,k)
load('XY.mat');
gx=X(find(lev==n));
gy=Y(find(lev==n));
m=sqrt(length(gx));
im=zeros(m,m);
for i=1:length(gx)
    im(gx(i),gy(i))=f(i+k);
end

imagesc(im);colorbar;