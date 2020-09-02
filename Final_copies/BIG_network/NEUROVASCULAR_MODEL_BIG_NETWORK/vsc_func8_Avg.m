function dPO2tisdt=vsc_func8_Avg(OEn,CMRO2n,L0,bndry,nE)
eA=bndry(5); eV=bndry(6);

dPO2tisdt=zeros(nE,1);
L=L0(eA+1:eV-1);
rho=30;
alpha=1.275e-12;
kk=alpha*pi*(rho^2)*L*1e-9;
dPO2tisdt(eA+1:eV-1)=(OEn-CMRO2n)./kk;

end