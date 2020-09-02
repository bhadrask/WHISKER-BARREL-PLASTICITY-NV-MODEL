function [Z] = get_constraint(A,lissom_size)
[x,y]=meshgrid(-20:20); xx=x(:);
x=0*x(:);y=y(:);
% figure; plot(x,y,'.')
a=2;
r = sqrt(a^2+x.^2+y.^2);
theta=atan(y./(x.^2+a.^2));
u=log(r);
v=atan(y);
% u=log(sec(v));; v=theta;
% v=
% A=2;
% figure; plot(A*u,A*v,'.r'); 
% hold on; plot(x,y,'.')
t=1:length(u);
f1=spline(t,A*u,1:0.05:length(u));
f2=spline(t,A*v,1:0.05:length(u));
 f2=f2- min(f2);
N = lissom_size(1);
xi = linspace(min(f1),max(f1),N) ;
yi = linspace(min(f2),max(f2),N) ;
[X,Y] = meshgrid(xi,yi) ;
%% get points lying inside polygon
idx = inpolygon(X(:),Y(:),f1,f2) ;
Z = zeros(size(X)) ;
Z(idx) = 1 ;  % assign one which lie inside 
%% plot to check
figure
imagesc(Z);
