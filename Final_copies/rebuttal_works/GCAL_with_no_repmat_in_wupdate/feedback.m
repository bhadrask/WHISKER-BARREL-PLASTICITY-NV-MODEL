function [lissom] = feedback(lissom,ii,opts)
rf1=lissom.layers{ii}.rf_fb(1)-1;
str=lissom.layers{ii+1}.stride;
cx=(repmat(1+(rf1)/2:str:lissom.layers{ii}.dim(1)-((rf1)/2),1,lissom.layers{ii+1}.dim(1)));
cz=(repmat(1+(rf1)/2:str:lissom.layers{ii}.dim(2)-((rf1)/2),lissom.layers{ii+1}.dim(2),1));
cy=reshape(cz,1,lissom.layers{ii+1}.dim(1)*lissom.layers{ii+1}.dim(2));
C=[cx;cy]';
G1=gaussian_filter_wt(lissom.layers{ii}.rf_fb(1),ceil(lissom.layers{ii}.rf_fb(1)/2));
cnt=0;
Mf=zeros(lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
Basemat=ones(lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
Z2d=reshape(lissom.layers{ii+1}.Zold,lissom.layers{ii+1}.dim(1),lissom.layers{ii+1}.dim(2));
for j=1: lissom.layers{ii+1}.dim(1)
    for i=1: lissom.layers{ii+1}.dim(2)
        cnt=cnt+1;
        G2=zeros(lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
      
        G2(cy-rf1/2:cy+rf1/2,cx-rf1/2:cx+rf1/2)=G1;
        G=G2(1:lissom.layers{ii}.dim(1),1:lissom.layers{ii}.dim(1)); 
        v=Z2d(cnt);
        cy=C(cnt,1);cx=C(cnt,2);
        X=v*Basemat;%repmat(v,lissom.layers{ii}.dim(1),lissom.layers{ii}.dim(2));
        Mf = Mf + opts.s*X.*G;
%         figure(1);imagesc(G); pause(0.3);
       
    end
end
lissom.layers{ii}.Zfb=reshape(Mf,lissom.layers{ii}.dim(1)*lissom.layers{ii}.dim(2),1);
%  figure(1);imagesc(Mf); pause(0.3);
end