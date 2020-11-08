 clc
clear;
close all;

load('output_4_plot_saturation_0.7_with_lesion.mat'); 
load('X_test.mat');
X1=X1/max(max(X1(:,:)));
dim1=lissom.layers{2}.dim;
[neu_size]=dim1;
neu_length=neu_size(1)*neu_size(2);

Xtest=X1;
ATP=x(end,3*nE+1:3*nE+neu_length);
min(ATP)


[yy,mean_act]=inst_response(lissom,opts,Xtest,label,ATP);
mean_act

% Xax=[0.2 0.3 0.31 0.32 0.33 0.35 0.45 0.5 1];
% 
% Yax=[-4 -4 1.14 1.28 1.28 1.28 1.28 1.28 1.28 ];
% plot(Xax,Yax);
xax=[0.2 0.5 0.55 0.6 0.62 0.63 0.65, 0.93];
yax=[1 1 1 1.13 1.29 1.67 1.67 1.67];
plot(xax,yax); xlabel('Saturation of oxygen at the inlet'); ylabel('D/C ratio');title('D/C ratio in the lesioned model with variation of inlet oxygen saturation')