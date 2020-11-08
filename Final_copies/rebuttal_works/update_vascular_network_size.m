clear;
close all;
clc;
Number_layers=11; %usually 7, cappilaries at the middle, 4th. 
mid=(Number_layers+1)/2;

split=4;
% split=4;

% layers=[1 16 256 4096 256 16 1];
 layers=[1 4 16 64 256 1024 256 64 16 4 1];
tot=sum(layers);
d0=zeros(tot,1);
dart(1)=8;
dvein(1)=8;
for vd=1:(numel(layers)/2)+1
    dart(vd+1)=dart(vd)+0.2*dart(vd);
    dvein(vd+1)=dvein(vd)+0.25*dvein(vd);
end
diam=[sort(dart(2:end),'descend'),dvein];
% diam=[30.5,19.5,12.5,8,15,23.4,36.6];

leng=[100*ones(1,length(dart(2:end))),250,100*ones(1,length(dart(2:end)))];
idx=1;

idx2=2;
idx3=1;
c3=1;
c=1;
N1=zeros(tot,1);
N2=zeros(tot,1);

for i=1:numel(layers)
    d0(idx:idx+layers(i)-1,1)=diam(i)*ones(layers(i),1);
        L0(idx:idx+layers(i)-1,1)=leng(i)*ones(layers(i),1);
 if i~=mid
    for j=1:layers(i)
      
          if  i<mid
        N1(idx2:idx2+split-1)=(c)*ones(split,1);
        N2(idx3)=c3;
        c=c+1;
        c3=c3+1;
         idx2=idx2+split;
           idx3=idx3+1;
          elseif i>mid
              N1(idx2)=c;
               N2(idx3:idx3+split-1)=(c3)*ones(split,1);
                 c3=c3+1;
              c=c+1;
              idx2=idx2+1;
              idx3=idx3+split;
          end
         
    end
 end
       
       
         
    
     idx=idx+layers(i);
end
N2(end)=N2(end-1)+1;
E=(1:tot)';
bndry=[layers(1)+layers(2)+layers(3)+layers(4)+1,layers(1)+layers(2)+layers(3)+layers(4)+layers(5)+1,layers(1)+layers(2)+layers(3)+layers(4)+layers(5)+1,layers(1)+layers(2)+layers(3)+layers(4)+layers(5)+layers(6),layers(1)+layers(2)+layers(3)+layers(4)+layers(5),layers(1)+layers(2)+layers(3)+layers(4)+layers(5)+layers(6)+1];
para=[60,25,10];
save('vessel_charac.mat','E','L0','N1','N2','bndry','d0','para','layers');
cn=1;
ct=1;
% for k=1:numel(layers);
%     dx=1:layers(k);
%     dim=sqrt(layers(k));
%     sq=reshape(dx,dim,dim);
%     for o=1:layers(k)
%         [Y(cn),X(cn)]=find(sq==o);
%         lev(cn)=ct;
%         cn=cn+1;
%     end
%    
%     ct=ct+1;
% end
% X=X';
% Y=Y';
% lev=lev';
% save('XY.mat','X','Y','lev');