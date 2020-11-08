clear;
close all;
clc;
load('lissom_10atpmax.mat');
load('amp5_X_test_sigma_4.mat');
load('Z_constraint');
X2=X1(:,:,12);
ATP=zeros(lissom.layers{2}.dim(1), lissom.layers{2}.dim(2));
inst_response(lissom,opts,X2,Z_constraint,label,ATP,ATPmax);
