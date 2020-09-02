clear all;
close all;clc;
i=0;
c=0;
count=0;
bb=6;
for ii=[6,18]
%     if ii==bb
%         yax=10:8:25;
%     else
        yax=[5:7:19];
%     end
    for  jj=yax
        c=c+1;
        center = [ii jj]; % again, example values
        %center=[0 0];
        for amp=1%:10:25
            count=count+10;
            sigma = [4 4];
            R = max(sigma(:));
            %[xgrid, ygrid] = meshgrid( center(1)-R: center(1)+R,center(2)-R: center(2)+R);
            [xgrid, ygrid] = meshgrid(1:24);
            for theta = 180%:-15:45
                
                %for theta=135:-90:45
                x = (((xgrid-center(1))*cosd(theta) + (ygrid-center(2))*sind(theta))/sigma(1)).^2;
                y = (((xgrid-center(1))*sind(theta) - (ygrid-center(2))*cosd(theta))/sigma(2)).^2;
                G = amp*exp(-(x+y)/4);i=i+1;
                 G(G<1e-10)=0;
                X1(:,:,i) = G;
                label(i)=count;
                % figure(1);
                % contour(xgrid, ygrid, G);pause(0.01);
                % figure(1);
                imagesc(G);
                max(max(G))
                pause(0.5)
            end
            
            [ii,jj]
        end
        
        
    end
    
end
  save('X_test.mat','X1','label');