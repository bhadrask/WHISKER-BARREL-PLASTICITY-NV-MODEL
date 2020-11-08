 function [res] = grating(stimSz, Orientation, Sf,  Tf, sinContrast, sinPhase)
    % Sf = Spatial frequency = no of sinusoidals per frame.
    % Tf = Temporal frequency = grating speed
    % sinContrast = more value increases the contrast
    y = [1:stimSz(1)] - (floor(stimSz(1)/2)+1);
    x = [1:stimSz(2)] - (floor(stimSz(2)/2)+1);
    t = [0:stimSz(3)-1];
    [y, x, t] = ndgrid(x, y, t);  % y = -y;
    Z =  sind(Orientation).*x + cosd(Orientation).*y;%rotate x and y values according to the desired orientation
%     imagesc(Z(:,:,1));
%     Z=0.5*ones([x y]); 
    res = sind(2*pi*Sf*Z -2*pi*Tf*t + 2*pi*sinPhase); %compute sin and cos functions
    res = sinContrast.*res./2 + .5;   %rescale to the range [0 1]
%     for lx=1:stimSz(1)
%         for ly=1:stimSz(2)
%             for lz=1:stimSz(3)
%                 if res(lx,ly,lz)< 2.4932
%                     res1(lx,ly,lz)=-50;
%                 end
%                 if res(lx,ly,lz)>= 2.4932
%                     res1(lx,ly,lz)= 50;
%                 end
%             end
%         end
%     end
 end
 
 % Call function as
 % x=grating([64  64  10],45,   5 , 5,1,0);