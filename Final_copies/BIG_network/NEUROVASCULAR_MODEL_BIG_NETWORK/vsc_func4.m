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
if numel(find(isnan(Fin))>0)>0
        disp('Fin gone wrong')
end
end

dVdtv(n1:eV-1)=Fin(n1:eV-1)-Fout(n1:eV-1);
 loc=find((imag(dVdtv))~=0);
 
if numel(loc)>0
       disp('dVdtv gone wrong')
end
