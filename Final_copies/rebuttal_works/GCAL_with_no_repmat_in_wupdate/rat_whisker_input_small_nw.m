clear all;
close all;clc;
i=0;
c=0;
count=0;
bb=6;
for ii=[29,56]
%     if ii==bb
%         yax=10:8:25;
%     else
        yax=[25:18:25+2*18];
%     end
    for  jj=yax
        c=c+1;
        center = [ii jj]; % again, example values
        %center=[0 0];
        count=count+10;
        for amp=5%30:2:34
            
            sigma = [10 10];
            R = max(sigma(:));
            %[xgrid, ygrid] = meshgrid( center(1)-R: center(1)+R,center(2)-R: center(2)+R);
            [xgrid, ygrid] = meshgrid(1:83);
            for theta = 180%:-15:45
                
                %for theta=135:-90:45
                x = (((xgrid-center(1))*cosd(theta) + (ygrid-center(2))*sind(theta))/sigma(1)).^2;
                y = (((xgrid-center(1))*sind(theta) - (ygrid-center(2))*cosd(theta))/sigma(2)).^2;
                G = amp*exp(-(x+y)/4);i=i+1;
                 G(G<1e-6)=0;
                X1(:,:,i) = G;
                label(i)=count;
                % figure(1);
                % contour(xgrid, ygrid, G);pause(0.01);
                % figure(1);
                 imagesc(G); h=colorbar;caxis([0 5]);pause(0.05);
                
            end
            
            [ii,jj]
        end
        
        
    end
    
end
  save('X_test_4vasc_amp5.mat','X1','label');