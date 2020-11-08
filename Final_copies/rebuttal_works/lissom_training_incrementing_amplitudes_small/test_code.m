clear; close all; clc;
[x,y]=meshgrid(-20:20); xx=x(:);
x=0*x(:);y=y(:);
figure; plot(x,y,'.')
a=1;
r = sqrt(a^2+x.^2+y.^2);
theta=atan(y./(x.^2+a.^2));
u=log(r);
v=atan(y);
% u=log(sec(v));; v=theta;
% v=
A=2;
figure; plot(A*u,A*v,'.r'); 
hold on; plot(x,y,'.')
t=1:length(u);
f1=spline(t,A*u,1:0.05:length(u));
f2=spline(t,A*v,1:0.05:length(u));
 f2=f2- min(f2);
figure();
plot(f1,f2)