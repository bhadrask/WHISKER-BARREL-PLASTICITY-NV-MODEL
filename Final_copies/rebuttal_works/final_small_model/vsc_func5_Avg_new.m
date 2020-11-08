function [dSdt,OE1,PO2_new]= vsc_func5_Avg_new(S,E,N1,N2,V,Fin,Fout,L0,d,PO2tisv,Sen)

kw=3.506*1e-18;
nE=length(E);
nN=N2(end)+1;
dSdt=zeros(nE,1);
A=zeros(nN-2,nN-2); b=zeros(nN-2,1); b1=zeros(nN-2,1);b2=zeros(nN-2,1);
Foutab=abs(Fout);
Fosm=zeros(nN-2,1);
% keyboard;
for i=1:nE-1
    if Fout(i)>0
        b1(N2(i))=b1(N2(i))+Foutab(i)*S(i);
    elseif Fout(i)<0
        Fosm(N2(i))=Fosm(N2(i))+Foutab(i);
    end
    ind=find(N1==N2(i));
    Finx=Fin(ind);
    nn=find(Finx<0);
    b2(N2(i))=sum(abs(Finx(Finx<0)).*S(ind(nn)));
    A(N2(i),N2(i))=sum(Finx(Finx>0))+Fosm(N2(i));
    b(N2(i))= b1(N2(i))+b2(N2(i));
end
% if det(A)==0
%     keyboard;
%     
% end
Sn=[Sen;A\b];

Sn(nN)=S(nE);
if numel(find(Sn>1))~=0 || numel(find(Sn<0))~=0
%     keyboard;
    disp('Sn gone wrong')
end
Sn(isnan(Sn))=0;
Sn(isinf(Sn))=0;
% Sn(Sn>1)=0;
% [cc1,cc2]=getGlobalx;
% if cc2==0
%     if numel(find(Sn>1))~=0 || numel(find(Sn<0))~=0
% %     keyboard;
%     disp('Sn gone wrong once')
%     end
%     cc2=cc2+1;
% end
% if cc2>=1
% if numel(find(Sn>1))~=0 || numel(find(Sn<0))~=0
% %     keyboard;
%     disp('Sn gone wrong again')
%    
%     cc2=cc2+1;
% 
% 
%  else
%      cc2=0;
% end
% if cc2>=7 || numel(find(isfinite(Sn)>0))~=numel(Sn)
%      disp(['Sn exploded  ',num2str( numel(Sn)-numel(find(isfinite(Sn)>0))), 'values went unbounded']);
%     error('STROKE OCCURED');
%  
% end
%     
%  end
% setGlobalx(cc1,cc2);

% if numel(find(Sn>1))~=0
%     keyboard;
% end

for i=1:nE
    PO2=exp(0.385*log((1/S(i)-1)^-1)+3.32-(72*S(i))^-1-(S(i)^6)/6);
%     if imag(PO2)~=0
%         keyboard;
%     end
    PO2_new(i)=PO2;
    
    OE=1*(kw*L0(i)*d(i))*(PO2-PO2tisv(i));
    
    if   OE<0
        OE=0;
        
    end
  
    OE11(i)=OE;
    loc1=find((imag(Fin))~=0);
    loc2=find((imag(Fout))~=0);
    if numel(loc1)>0 ||numel(loc2)>0
        keyboard;
    end
    if Fin(i)>=0
        Sin=Sn(N1(i)+1);
    elseif Fin(i)<0
        Sin=S(i);
    end
 
    if Fout(i)>=0
        Sout=S(i);
    elseif Fout(i)<0
        Sout=Sn(N2(i)+1);
    end
    
    dSdt(i)=(Fin(i)*Sin-Fout(i)*Sout-OE*((4*2.3*1e-9)^-1)-(Fin(i)-Fout(i))*S(i))/V(i);
    
    
    
end
dSdt(isnan(dSdt))=0;
dSdt(isinf(dSdt))=0;
OE1=OE11';
