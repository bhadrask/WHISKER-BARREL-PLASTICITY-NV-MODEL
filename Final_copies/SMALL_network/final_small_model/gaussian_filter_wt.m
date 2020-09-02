function [W]=gaussian_filter_wt(k,sd)

[ygrid, xgrid] = meshgrid(1:k);
ci=(size(xgrid,1)+1)/2;
cj=ci;

        x = (((xgrid-ci) + (ygrid-cj))/sd).^2;
        y = (((xgrid-ci) - (ygrid-cj))/sd).^2;
        w= exp(-(x+y)/4);
        w(w<1e-4)=0;
        W=normalize(w,1);
        