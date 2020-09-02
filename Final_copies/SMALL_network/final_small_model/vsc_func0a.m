function [Ae,Pn,Pe,Fin,Fout]=vsc_func0a(E,N1,N2,bndry,R,V,para)

nE=length(E);
nN=N2(end)+1;
Pa=para(1);Pv=para(2);Pic=para(3);
Pe=zeros(nE,1);
Fin=zeros(nE,1);
Fout=zeros(nE,1);
A=zeros(nN-2);
b=zeros(nN-2,1);

for i=1:nE
    if N1(i)==0
        A(N2(i),N2(i))=A(N2(i),N2(i))-1/R(i);
        b(N2(i))=b(N2(i))-Pa/R(i);
    elseif N2(i)==N2(end)
        A(N1(i),N1(i))=A(N1(i),N1(i))-1/R(i);
        b(N1(i))=b(N1(i))-Pv/R(i);        
    else
        A(N1(i),N1(i))=A(N1(i),N1(i))-1/R(i);
        A(N2(i),N2(i))=A(N2(i),N2(i))-1/R(i);
        A(N1(i),N2(i))=A(N1(i),N2(i))+1/R(i);
        A(N2(i),N1(i))=A(N2(i),N1(i))+1/R(i);
    end
end
Pn=[Pa;(A\b);Pv];
for i=1:nE
%     i
%    n1= (N1(i)+1)
%     n2=(N2(i)+1)
    Pe(i)=mean([Pn(N1(i)+1),Pn(N2(i)+1)]);
    Fin(i)=(Pn(N1(i)+1)-Pe(i))/(R(i));
    Fout(i)=(Pe(i)-Pn(N2(i)+1))/(R(i));
% Fout(i)=Fin(i);
%     pe=Pe(i)
%     fi=Fin(i)
%     fo=Fout(i)
   del(i)=(Fin(i)-Fout(i));
end

eA=bndry(5);
% beta=2*ones(nE-eA,1); 
% Ae=V(eA+1:end)./((Pe(eA+1:end)-Pic).^(1./beta));

%  Ae=Pe(eA+1:end)-Pic;
 Ae=(Pe(1:end)-Pic)./((V(1:end)).^1.7);
% Ar=(Pe(1:eA)-Pic)./V(1:eA);
end