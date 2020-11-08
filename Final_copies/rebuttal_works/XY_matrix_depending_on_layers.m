% load('XY.mat');
clear;
Number_layers=11; %usually 7, cappilaries at the middle, 4th. 
mid=(Number_layers+1)/2;

split=4;
% split=4;

% layers=[1 16 256 4096 256 16 1];
 layers=[1 4 16 64 256 1024 256 64 16 4 1]; %usually 7, cappilaries at the middle, 4th.
tot=sum(layers);

mid=(Number_layers+1)/2;
split2d=sqrt(split);
m=1;
X1(m)=1;Y1(m)=1;
for i=1:mid-1
    ele=layers(i)*split;
    basic=zeros(sqrt(ele),sqrt(ele));
    c=1;d=1;
    for j=1:layers(i)
        basic=zeros(sqrt(ele),sqrt(ele));
        basic(c:c+split2d-1,d:d+split2d-1)=1;
        figure(1);
        imagesc(basic);pause(0.2);
        basic;
        [locx,locy]=find(basic==1);
        X1(m+1:m+split)=locy;
        Y1(m+1:m+split)=locx;
        m=m+split;
        
        if rem(j,split2d)==0
            if rem(j,split)~=0
                d=d+split2d;
                c=c-split2d*(split2d-1);
                
            elseif rem(j,split)==0 && rem(j,sqrt(ele))~=0
                d=d-split2d*(split2d-1);
                c=c+split2d;
            elseif rem(j,split)==0 && rem(j,sqrt(ele))==0
                c=1;
                d=d+split2d;
            end
        else
            c=c+split2d;
        end
        %         R=X1';
        %    S=X(1:length(R));
   
    end
    
end
n=1;
bal=tot-length(X1);
X3={};
Y3={};
cn=length(numel(layers):-1:mid+1);
for k=numel(layers):-1:mid+1
    X2=X1(n:n+layers(k)-1);
    Y2=Y1(n:n+layers(k)-1);
    X3{cn}=X2;
    Y3{cn}=Y2;
    n=n+layers(k);
    cn=cn-1;
end
for p=1:numel(X3)
    X1=[X1 X3{p}];
    Y1=[Y1 Y3{p}];
end
cn=1;
ct=1;
for k=1:numel(layers);
    
    for o=1:layers(k)
        
        lev(cn)=ct;
        cn=cn+1;
    end
    
    ct=ct+1;
end

X=X1';
Y=Y1';
lev=lev';
save('XY.mat','X','Y','lev');
find(X1>64)