clear all;
close all;clc;
i=0;
c=0;
count=0;
bb=18;
for ii=bb:10:60
    if ii==bb
        yax=20:10:55;
    else
        yax=14:10:60;
    end
    for  jj=yax
        c=c+1;
        center = [ii jj]; % again, example values
        %center=[0 0];
        for amp=15;%5:10:40
            count=count+1;
            sigma = [6 6];
            R = max(sigma(:));
            %[xgrid, ygrid] = meshgrid( center(1)-R: center(1)+R,center(2)-R: center(2)+R);
            [xgrid, ygrid] = meshgrid(1:70);
            for theta = 180%:-15:45
                
                %for theta=135:-90:45
                x = (((xgrid-center(1))*cosd(theta) + (ygrid-center(2))*sind(theta))/sigma(1)).^2;
                y = (((xgrid-center(1))*sind(theta) - (ygrid-center(2))*cosd(theta))/sigma(2)).^2;
                G = amp*exp(-(x+y)/2);i=i+1;
                 G(G<1e-6)=0;
                X1(:,:,i) = G;
                label(i)=count+10;
                % figure(1);
                % contour(xgrid, ygrid, G);pause(0.01);
                % figure(1);
                imagesc(G); h=colorbar;caxis([0 25]);pause(0.05);
                max(max(G))
             
            end
            
            [ii,jj]
        end
        
        
    end
    
end
 save('Final_X_test_sigma_6_AMP15.mat','X1','label');