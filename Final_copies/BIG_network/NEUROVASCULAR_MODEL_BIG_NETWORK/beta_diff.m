function [dbetadt]= beta_diff(neu, bndry,Ae,beta,W,filename,loc,N1,N2)

load(filename);
G=[X(find(lev==4)),Y(find((lev==4)))];
eA=bndry(5); eV=bndry(6);
caps=sqrt(eV-eA-1);

neu_net=zeros(1,size(G,1));
neu_net1=conv2(neu,W,'same');
neu_net2=reshape(neu_net1,caps*caps,1);
neu_net(loc)=neu_net2;

neu_net=neu_net';
Ae0=Ae(1:eA);
Ae1=Ae(eA+1:eA+length(neu_net),1);
beta_t=beta(eA+1:eA+length(neu_net),1);
Ae2=Ae(eA+length(neu_net)+1:end,1);
dbetadt1=(((-beta_t)+Ae1.*(1-(neu_net/10))))/0.5; % c1,t1

B=zeros(eA,1);

for i=eA+1:eV-1
    B(N1(i))=B(N1(i))+neu_net(i-eA);
end
B=B/numel(find(N1==N2(7)));
B0=beta(1:eA);
dB0dt=(((-B0)+Ae0.*(1-(B)/1)))/1;
dbetadt=[dB0dt;dbetadt1;zeros(size(Ae2))];