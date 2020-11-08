function dxdt=function4DDE(t,x,z,Ae,opts,input1,CMRO2nref,ATPref,PO2tis4art,W_gauss,loc_vas,sen,file_xy,file_vess_char,training,LISSOM,tau,Fin0,Fout0,V0,PO2tismax,ATPmax,thr_atp)
%% definitions

load(file_xy);
load(file_vess_char);
%  [lissom]=getGlobalx;
nE=length(E);
neu_size=LISSOM.layers{2}.dim;
neu_length=neu_size(1)*neu_size(2);
% ATP_ch=x(3*nE+1:3*nE+neu_length);
% ATP_check=min(ATP_ch)
% if  t>0.51
%     keyboard;
% end
tstart=0.5;
tstop=1.5;
 if  t>tstart && t<=tstop
     input=input1;%keyboard;
 else
  
     input=zeros(size(input1,1),size(input1,2));
%      keyboard;
    
 end
 t1=t-tau;
  if  t1>tstart && t1<=tstop
%        keyboard;
     input2=input1;
 else
  
     input2=zeros(size(input1,1),size(input1,2));
%      keyboard;
%     
  end


%% Loading previous iteration values to find the increments
V=x(1:nE);
S=x(nE+1:2*nE);
PO2tis=x(2*nE+1:3*nE);
ATP=x(3*nE+1:3*nE+neu_length);
% if min(ATP)<2.1
%     keyboard;
% end
ATP_prev=z(3*nE+1:3*nE+neu_length);
beta=x(3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n=x(3*nE+neu_length+nE+1:end);
% neu_prev=getGlobalneu;
%% Adjusting the values if something has gone beyond its given range
locv=find(V<0);
V(locv)=0;
x(locv)=0;
locsat=find(S<0);
S(locsat)=0;
x(nE+locsat)=0;
locsp=find(PO2tis<0);
PO2tis(locsp)=0;
x(2*nE+locsp)=0;
locs=find(ATP<0);
ATP(locs)=0;
x(3*nE+locs)=0;

 ATP_prev(ATP_prev<0)=0;


 

%% Calculations to find updated dependent variables
newd=2*sqrt((V*1e9)./(pi*L0));
R=128*15*(L0/2)*1e3./(pi*newd.^4);
PO2tisv=project2vasc(bndry,PO2tis,PO2tis4art,W_gauss,loc_vas);

 lissom=getGlobaly;
lissom1=lissom;
lissomref=LISSOM;
Zold_ref=lissomref.layers{2}.Zold;
[lissom.layers{2}.Zold,lissom1.layers{2}.Zold]=getGlobalx;
 [cnt,neuout,ntime]=getGlobalneu;
%% Neural and vascular variables increment
% neu=inst_response(lissom,opts,input,ATP,ATPref);

if (ntime(cnt)<=tstart && t>tstart)||(ntime(cnt)<tstop && t>=tstop)
    compulsory=1;
else
    compulsory=0;
end
if ntime(cnt)>t
    [~,idx]=min(abs(ntime-t));
    cnt=max(idx);
     lissom.layers{2}.Zold=neuout(:,cnt);
 end
 if (abs(ntime(cnt)-t)<0.1 && compulsory==0) %t2
     lissom.layers{2}.Zold=neuout(:,cnt);
     neu=lissom.layers{2}.Zold;
     
 else
      kcnt=round(10*abs((abs(ntime(cnt)-t))));
    if kcnt==0
        keyboard;
    end
    
    for i=1:kcnt
      
       if numel(find(input>0))>0
       [neu,lissom]=inst_response_updat_depending_vsclr_fdk(lissom,opts,input,ATP,training,ATPmax,thr_atp);
       setGlobaly(lissom);
       lissom.layers{2}.Zold=neu;
       else
             lissom.layers{2}.Zold=neuout(:,cnt);
     neu=lissom.layers{2}.Zold;
       end
%        keyboard;
%         mno=cell2mat(lissom.layers{2}.waff);
%         diff_max=max(max(mno-Waffref))
%         diff_min=min(min(mno-Waffref))
    end
%     disp(['kcnt =',num2str(kcnt),'  maxneu=',num2str(max(max(neu))),'  time = ',num2str(t)] );

    cnt=cnt+1;
    neuout(:,cnt)=neu;
    ntime(cnt)=t;
    setGlobalneu(cnt,neuout,ntime);

 end

 
if t1>0
   [~,idx2]=min(abs(ntime-t1));
    cnt2=max(idx2);
 lissom1.layers{2}.Zold=neuout(:,cnt2);
    neu_prev=lissom1.layers{2}.Zold;

 else
     lissom1.layers{2}.Zold=Zold_ref;
     neu_prev=lissom1.layers{2}.Zold;
 end

 setGlobalx(neu,neu_prev);
dbetadt=beta_diff(neu_prev, bndry,Ae,beta,W_gauss,file_xy,loc_vas,N1,N2,layers);
% [Pnv,Pe,dVdt]=vsc_func2(E,N1,N2,bndry,R,V,para,beta);
% [PnA,Pe]=vsc_func3(E,N1,N2,bndry,R,Pe,dVdt,para);
[Pnv,Pe,dVdt]=vsc_func2_new(E,N1,N2,bndry,R,V,para,beta,lev,layers);
[PnA,Pe]=vsc_func3_new(E,N1,N2,bndry,R,Pe,dVdt,para,lev,layers);
Pn=[PnA;Pnv];
Pnold=Pn;
% if numel(find(Pn>60))~=0
%     keyboard;
% end
[dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe,Fin0,Fout0,V0);
[dSdt,OE1,PO2_new]= vsc_func5_Avg_new(S,E,N1,N2,V,Fin,Fout,L0,newd,PO2tisv,sen);
[dCMRO2ndt,OEn]=vsc_func7_Avg_bsk_mod(OE1,bndry,neu,CMRO2nref,W_gauss,file_xy,CMRO2n,PO2tis,PO2tismax,layers);
% [dATPdt,Jatp]=vsc_func9_Avg(CMRO2n,ATP,neu,lissom,CMRO2nref,ATPref,OE0,OEn);
[dATPdt,Jatp]=vsc_func9_Avg_new(CMRO2n,ATP,neu,lissom,CMRO2nref,ATPref,PO2tis,PO2tismax,bndry);
dPO2tisdt=vsc_func8_Avg(OEn,CMRO2n,L0,bndry,nE);
dVdt=dVdtv;
dxdt=[dVdt;dSdt;dPO2tisdt;dATPdt;dbetadt;dCMRO2ndt];
if numel(find(isnan(dxdt)))>0
    keyboard;
end
disp(t)
