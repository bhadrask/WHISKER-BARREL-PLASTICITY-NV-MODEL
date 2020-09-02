function [dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe,Fin0,Fout0,V0)
% function [dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe)
eA=bndry(5);eV=bndry(6);
n1=bndry(1);
nE=length(E);
Fin=zeros(nE,1);
Fout=zeros(nE,1);
dVdtv=zeros(nE,1);
for i=1:nE
    Fin(i)=(Pn(N1(i)+1)-Pe(i))/R(i);
    Fout(i)=(Pe(i)-Pn(N2(i)+1))/R(i);
%     if (Fin(i)<0)||(Fout(i)<0)
%         keyboard;
%     end
% if abs(Fin(i))<1e-10
%     Fin(i)=0;
% end
% if abs(Fout(i))<1e-10
%     Fout(i)=0;
% end
% 
% if numel(find(isnan(Fin))>0)>0
%     keyboard;
% end
end
Fin(isnan(Fin))=0;
Fin(isinf(Fin))=0;
Fout(isnan(Fout))=0;
Fout(isinf(Fout))=0;
const=1e-4;
dVdtv(n1:eV-1)=Fin(n1:eV-1)-Fout(n1:eV-1);
% dVdtv(eA+1:eV-1)=V0(eA+1:eV-1).*(sign(Fin(eA+1:eV-1)-Fout(eA+1:eV-1)).*abs((Fin(eA+1:eV-1)-Fout(eA+1:eV-1))).^0.4)./(abs(Fin0(eA+1:eV-1))).^0.4;
% dVdtv(eA+1:eV-1)=const*V0(eA+1:eV-1).*((Fin(eA+1:eV-1)-Fout(eA+1:eV-1))./(Fin0(eA+1:eV-1)-Fout0(eA+1:eV-1))).^0.4;
% dVdtv(eA+1:eV-1)=const*V0(eA+1:eV-1).*sign((Fin(eA+1:eV-1)./Fin0(eA+1:eV-1))).*(abs((Fin(eA+1:eV-1)./Fin0(eA+1:eV-1))-(Fout(eA+1:eV-1)./Fout0(eA+1:eV-1)))).^0.4;
%   keyboard
dVdtv(1);
% ll=find(abs(dVdtv)==max(abs(dVdtv)));
% vm=max(abs(dVdtv))
% fi=Fin(eA+ll)
% fo=Fout(eA+ll)
 loc=find((imag(dVdtv))~=0);
if numel(loc)>0
    keyboard
Pe
R
    error('non invertible in vsc_func4');
end
