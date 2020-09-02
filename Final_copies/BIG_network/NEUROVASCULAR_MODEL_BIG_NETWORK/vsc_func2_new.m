function [Pnv,Pe,dVdt]=vsc_func2_new(E,N1,N2,bndry,R,V,para,beta)
n1=bndry(1);
n2=bndry(2);
nE=length(E);
nN=N2(end)+1;
dVdtv=zeros(nE,1);
Pv=para(2);Pic=para(3);
A=zeros(nN-1-n1,nN-1-n1); b=zeros(nN-1-n1,1); b1=zeros(nN-1-n1,1); b2=zeros(nN-1-n1,1);

Pe=zeros(nE,1);

Pe(n1:end)=((V(n1:end).^1).*(beta(n1:end)))+Pic;

for i=n1:nE
    if N1(i)>=2 && N1(i)<n1
        A(N2(i)-n1+1,N2(i)-n1+1)=A(N2(i)-n1+1,N2(i)-n1+1)+1/R(i);
        b(N2(i)-n1+1)=b(N2(i)-n1+1)+Pe(i)/R(i);
        
    elseif N1(i)>=n1 && N1(i)<n2 %N2(i)~=N2(end)
 
        A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
        b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
        temp=Pe(N1==N2(i))/R(N1==N2(i));
        temp2=1/R(N1==N2(i));
        b1(N2(i)-n1+1)= b1(N2(i)-n1+1)+Pe(i)/R(i);
        b(N2(i)-n1+1)=b1(N2(i)-n1+1)-temp;
        A(N2(i)-n1+1,N1(i)-n1+1)= A(N2(i)-n1+1,N1(i)-n1+1)+1/R(i);
        A(N2(i)-n1+1,N2(i)-n1+1)=-temp2;
    elseif N1(i)>=n2 && N2(i)~=N2(end)
        A(N2(i)-n1+1,N1(i)-n1+1)=A(N2(i)-n1+1,N1(i)-n1+1)+1/R(i);
        temp4=1/R(N1==N2(i));
        A(N2(i)-n1+1,N2(i)-n1+1)=-temp4;
        temp3=Pe(N1==N2(i))/R(N1==N2(i));
        b2(N2(i)-n1+1)= b2(N2(i)-n1+1)+Pe(i)/R(i);
        b(N2(i)-n1+1)= b2(N2(i)-n1+1)-temp3;
        
    end
end

loc=find((imag(A))~=0);
if numel(loc)>0
    keyboard;
    error('non invertible in vsc_func2');
end
Pnv=[(A\b);Pv];

dVdtv(n2:end)=0;
dVdt=dVdtv;
end

