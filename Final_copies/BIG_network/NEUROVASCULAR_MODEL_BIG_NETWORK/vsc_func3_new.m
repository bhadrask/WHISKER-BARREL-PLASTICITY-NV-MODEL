function [Pna,Pe]=vsc_func3_new(E,N1,N2,bndry,R,Pe,dVdt,para)

n1=bndry(1);
n2=bndry(2);

Pa=para(1);
A=zeros(n1-1);
b=zeros(n1-1,1);
keyboard;
for i=1:n2-1
    if N1(i)==0
        A(N2(i),N2(i))=A(N2(i),N2(i))-1/R(i);
        b(N2(i))=b(N2(i))+dVdt(i)-Pa/R(i);
    elseif (N1(i)>=2)&&(N1(i)<n1)
        A(N1(i),N1(i))=A(N1(i),N1(i))-2/R(i);
        b(N1(i))=b(N1(i))-2*Pe(i)/R(i);
        if N2(i)<n1
            A(N2(i),N2(i))=A(N2(i),N2(i))-2/R(i);
            b(N2(i))=b(N2(i))-2*Pe(i)/R(i);
        end
    else
        A(N1(i),N1(i))=A(N1(i),N1(i))-1/R(i);
    
        A(N2(i),N2(i))=A(N2(i),N2(i))-1/R(i);
        A(N1(i),N2(i))=A(N1(i),N2(i))+1/R(i);
        A(N2(i),N1(i))=A(N2(i),N1(i))+1/R(i);
        b(N1(i))=b(N1(i))+dVdt(i);
        b(N2(i))=b(N2(i))+dVdt(i);       
    end
end


Pna=[Pa;(A\b)];

for i=1:n1-1
    Pe(i)=(1/2)*(Pna(N1(i)+1)+Pna(N2(i)+1)-R(i)*dVdt(i));
end