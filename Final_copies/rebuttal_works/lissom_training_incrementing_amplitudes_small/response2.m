function [q]=response2(lissom,opts,X1,label)
% load trained network and the testdata with label if you don't want
%to run it as function
yy=[];zz=[];

a=X1 ;xdim = size(a(:,:,1));


% sel=steps;%1:steps:size(a,3);
% 
% p=sel;

p=1;
%for i=1:length(steps)
   
   % for inp=1:steps(i)
       % for loop=1:4
       for i1= 1: size(X1,3) 
            inpsel=a(:,:,i1); 
            figure(1);subplot(3,3,1);
            imagesc(inpsel);pause(0.1);
            for jj=2:numel(lissom.layers)
                [lissom] = rf_m(lissom,double(inpsel),jj,opts);
                [lissom] = activate(lissom,jj,opts);
                if jj>1 && jj<numel(lissom.layers) && opts.control==1
                    [lissom] = feedback(lissom,jj,opts);
                end
                [lissom] = lat_dynamics(lissom,jj,opts);
                if jj==2 
                    yy=cat(3,yy,lissom.layers{jj}.Zold);
                elseif jj==3 
                    zz=cat(3,zz,lissom.layers{jj}.Zold);
                end
                
                inpsel=lissom.layers{jj}.Zold;
                lissom.layers{jj}.Zabove=lissom.layers{jj}.Zold;
                
            end
           
        %end
        
        
        
   
     
    for oo=2:numel(lissom.layers)
        lissom.layers{oo}.Zold=zeros(lissom.layers{oo}.dim);
        lissom.layers{oo}.Zabove=zeros(lissom.layers{oo}.dim);
        lissom.layers{oo}.Z = zeros(lissom.layers{oo}.dim);
    end

    
    
end
    [gg,q]=max(yy,[],3)
   % figure; imagesc(labeldir(q));title('layer2 direction label');colorbar;
    figure; imagesc(label(q));title('layer2 orientation label');colorbar;
    [pp,loc]=max(zz,[],3);
   % figure;imagesc(labeldir(loc));title('layer3 direction label');colorbar;
     figure;imagesc(label(loc));title('layer3 orientation label');colorbar;
    %  loc=find(gg==0);
    %  q(loc)=length(label);
    %   figure; imagesc(label(q));colorbar;
    %   gg
    %   q
    %   I=label(q);
    %   Iblur2 = imgaussfilt(I,4);
    %   figure;imagesc(Iblur2)
end

% 


