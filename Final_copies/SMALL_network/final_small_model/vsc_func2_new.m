function [Pnv,Pe,dVdt]=vsc_func2_new(E,N1,N2,bndry,R,V,para,beta)
n1=bndry(1);
n2=bndry(2);
e1=bndry(3);
eA=bndry(5);
nE=length(E);
nN=N2(end)+1;
dVdtv=zeros(nE,1);
Pv=para(2);Pic=para(3);
A=zeros(nN-1-n1,nN-1-n1); b=zeros(nN-1-n1,1); b1=zeros(nN-1-n1,1); b2=zeros(nN-1-n1,1);
A1=zeros(nN-1-n1,nN-1-n1); 
Pe=zeros(nE,1);
% Peold=Pr(1:nE);
% Pnold=Pr(nE+1:end);
%  Pe(eA+1:end)=((V(eA+1:end)./V0(eA+1:end)).^beta).*Ae+Pic;
Pe(n1:end)=((V(n1:end).^1.7).*(beta(n1:end)))+Pic;
% Pe(25)
% beta(25)
%   keyboard;
% 
% for i=n1:nE
%     if N1(i)<n2
%         A(N2(i)-n1+1,N2(i)-n1+1)=A(N2(i)-n1+1,N2(i)-n1+1)+1/R(i);
%         b(N2(i)-n1+1)=b(N2(i)-n1+1)+Pe(i)/R(i);
% %         if N1(i)>=2 && N1(i)<n1
% %              A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
% %         b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
% %         A(N2(i)-n1+1,N2(i)-n1+1)=A(N2(i)-n1+1,N2(i)-n1+1)+1/R(i);
% %         b(N2(i)-n1+1)=b(N2(i)-n1+1)+Pe(i)/R(i);
% %         end
%     elseif N2(i)==N2(end)
%         A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
%         b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
%     else
%         A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
%         b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
%         A(N2(i)-n1+1,N2(i)-n1+1)=A(N2(i)-n1+1,N2(i)-n1+1)+1/R(i);
%         b(N2(i)-n1+1)=b(N2(i)-n1+1)+Pe(i)/R(i);
%     end
% end
for i=n1:nE
    if N1(i)>=2 && N1(i)<n1
        A(N2(i)-n1+1,N2(i)-n1+1)=A(N2(i)-n1+1,N2(i)-n1+1)+1/R(i);
        b(N2(i)-n1+1)=b(N2(i)-n1+1)+Pe(i)/R(i);
        
    elseif N1(i)>=n1 && N1(i)<n2 %N2(i)~=N2(end)
         A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
        b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
         temp=Pe(find(N1==N2(i)))/R(find(N1==N2(i)));
         temp2=1/R(find(N1==N2(i)));
        b1(N2(i)-n1+1)= b1(N2(i)-n1+1)+Pe(i)/R(i);
        b(N2(i)-n1+1)=b1(N2(i)-n1+1)-temp;
        A(N2(i)-n1+1,N1(i)-n1+1)= A(N2(i)-n1+1,N1(i)-n1+1)+1/R(i);
         A(N2(i)-n1+1,N2(i)-n1+1)=-temp2;
%     elseif N1(i)>=n2 && N2(i)~=N2(end)
%        
%         A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
%         b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
%         A(N2(i)-n1+1,N2(i)-n1+1)=A(N2(i)-n1+1,N2(i)-n1+1)+1/R(i);
%         b(N2(i)-n1+1)=b(N2(i)-n1+1)+Pe(i)/R(i);
   elseif N1(i)>=n2 && N2(i)~=N2(end)
        A(N2(i)-n1+1,N1(i)-n1+1)=A(N2(i)-n1+1,N1(i)-n1+1)+1/R(i);
          temp4=1/R(find(N1==N2(i)));
         A(N2(i)-n1+1,N2(i)-n1+1)=-temp4;
      
         temp3=Pe(find(N1==N2(i)))/R(find(N1==N2(i)));
        
         b2(N2(i)-n1+1)= b2(N2(i)-n1+1)+Pe(i)/R(i);
         b(N2(i)-n1+1)= b2(N2(i)-n1+1)-temp3;
%         elseif N2(i)==N2(end)
%         A(N1(i)-n1+1,N1(i)-n1+1)=A(N1(i)-n1+1,N1(i)-n1+1)+1/R(i);
%         b(N1(i)-n1+1)=b(N1(i)-n1+1)+Pe(i)/R(i);
    end
end

% det(A)
loc=find((imag(A))~=0);
if numel(loc)>0
    
  R
    
  V
  keyboard;
    error('non invertible in vsc_func2');
end
Pnv=[(A\b);Pv];
Pnv(isnan(Pnv))=0;
Pnv(isinf(Pnv))=0;
% Pn=[Pnold(1:(nN-length(Pnv)));Pnv];
% Pn=Pnold;
% Peo=[Peold(1:n2-1);Pe(eA+1:end)];
%   for i=1:nE
%     Fin(i)=(Pn(N1(i)+1)-Peo(i))/R(i);
%     Fout(i)=(Peo(i)-Pn(N2(i)+1))/R(i);
%     if (Fin(i)<0)||(Fout(i)<0)
%         keyboard;
%     end
%     
% end
% dVdtv(n2:end)=Fin(n2:end)-Fout(n2:end);
if numel(find(isnan(Pnv))>0)>0
    keyboard;
end
dVdtv(n2:end)=0;
dVdt=dVdtv;   
end
    
