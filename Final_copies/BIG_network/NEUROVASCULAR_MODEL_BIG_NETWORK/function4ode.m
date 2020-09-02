function dxdt=function4ode(t,x,Ae,lissom,opts,Amax,IN,CMRO2nref,ATPref,Z_constraint,PO2tis4art,W_gauss,loc_vas,sen,file_xy,file_vess_char,PO2tismax)
%% definitions
load(file_xy);
load(file_vess_char);
nE=length(E);
neu_size=lissom.layers{2}.dim;
neu_length=neu_size(1)*neu_size(2);
input1=IN.input;
input=input1;
tstart=IN.Tstart;
tstop=IN.Tstop;
tau=IN.Tau;
if (t<=tstart) || (t>tstop)
    input=zeros(size(input1));
end
t1=t-tau;
if (t1<=tstart) || (t1>tstop)
    input2=zeros(size(input1));
else
    input2=input1;
end
% if ((t>tstart)&& (t<tstart+0.001)) ||(t1>tstart)||(t>tstop)||(t1>tstop)
%     keyboard;
% end


%% Loading previous iteration values to find the increments
V=x(1:nE);
S=x(nE+1:2*nE);
PO2tis=x(2*nE+1:3*nE);
ATP=x(3*nE+1:3*nE+neu_length);
% ATP_prev=z(3*nE+1:3*nE+neu_length);
beta=x(3*nE+neu_length+1:3*nE+neu_length+nE);
CMRO2n=x(3*nE+neu_length+nE+1:end);
%% Adjusting the values if something has gone beyond its given range
% locsat=find(S<0);
% S(locsat)=0;
% x(nE+locsat)=0;
% locsp=find(PO2tis<0);
% PO2tis(locsp)=0;
% x(2*nE+locsp)=0;
% locs=find(ATP<0);
% ATP(locs)=0;
% x(3*nE+locs)=0;
%% Calculations to find updated dependent variables
newd=2*sqrt((V*1e9)./(pi*L0));
R=128*15*(L0/2)*1e3./(pi*newd.^4);
PO2tisv=project2vasc(bndry,PO2tis,PO2tis4art,W_gauss,loc_vas);
%% displaying non adjustable errors due to integration time step increase
if numel(find(S>1))~=0 || numel(find(S<0))~=0
    disp('warning S has gone wrong');
end

% keyboard;
%% Neural and vascular variables increment

lissom1=lissom;
lissomref=lissom;
Zold_ref=lissomref.layers{2}.Zold;
[lissom.layers{2}.Zold,lissom1.layers{2}.Zold]=getGlobalx;
  [cnt,neuout,ntime]=getGlobalneu;
 %% steps to integrate the ode with neural network time
%to check if input is in transition
if (ntime(cnt)<=tstart && t>tstart)||(ntime(cnt)<tstop && t>=tstop)
    compulsory=1;
else
    compulsory=0;
end
if ntime(cnt)>t
    [~,idx]=min(abs(ntime-t));
    cnt=max(idx);
     lissom.layers{2}.Zold=neuout(:,:,cnt);
 end
 if abs(ntime(cnt)-t)<0.1 && compulsory==0 %t2
     lissom.layers{2}.Zold=neuout(:,:,cnt);
     neu=lissom.layers{2}.Zold;
     
 else
      kcnt=round(10*abs((abs(ntime(cnt)-t))));
    if kcnt==0
        keyboard;
    end
    
    for i=1:kcnt
        neu=inst_response(lissom,opts,input,Z_constraint,ATP,Amax);
        lissom.layers{2}.Zold=neu;
    end
%     disp(['kcnt =',num2str(kcnt),'  maxneu=',num2str(max(max(neu))),'  time = ',num2str(t)] );

    cnt=cnt+1;
    neuout(:,:,cnt)=neu;
    ntime(cnt)=t;
    setGlobalneu(cnt,neuout,ntime);
%      neu=inst_response(lissom,opts,input,Z_constraint,ATP,Amax);
%       cnt=cnt+1;
%     neuout(:,:,cnt)=neu;
%     ntime(cnt)=t;
%      setGlobalneu(cnt,neuout,ntime);
 end

 if t1>0
   [~,idx2]=min(abs(ntime-t1));
    cnt2=max(idx2);
%     compulsory2=0;
%     if cnt2>1
%     if (ntime(cnt2-1)<=tstart && t1>tstart)||(ntime(cnt2-1)<tstop && t1>=tstop)
%         compulsory2=1;   
%     end
%     end
%     if abs(ntime(cnt2)-t1)<0.005 && compulsory2==0 %t2
    lissom1.layers{2}.Zold=neuout(:,:,cnt2);
    neu_prev=lissom1.layers{2}.Zold;
%     else
%         neu_prev=inst_response(lissom1,opts,input2,Z_constraint,ATP,Amax);
%     end
 else
     lissom1.layers{2}.Zold=Zold_ref;
     neu_prev=lissom1.layers{2}.Zold;
 end
% if t1<=tstart
%     lissom1.layers{2}.Zold=Zold_ref;
% end
% neu=inst_response(lissom,opts,input,Z_constraint,ATP,Amax);
% neu_prev=inst_response(lissom1,opts,input2,Z_constraint,ATP,Amax);
% if t1>tstart
  setGlobalx(neu,neu_prev);
% else
%     setGlobalx(neu,lissom1.layers{2}.Zold);
% end

% % if ismember(round(t,3),round(0:0.1:4,3))
% %     [cnt,neuout]=getGlobalneu;
% %     cnt=cnt+1;
% %     neuout(:,:,cnt)=neu;
% %      setGlobalneu(cnt,neuout);
% % end
%  cnt=cnt+1;
%     neuout(:,:,cnt)=neu;
%     ntime(cnt)=t;
%      setGlobalneu(cnt,neuout,ntime);
dbetadt=beta_diff(neu_prev, bndry,Ae,beta,W_gauss,file_xy,loc_vas,N1,N2);
[Pnv,Pe,dVdt]=vsc_func2_new(E,N1,N2,bndry,R,V,para,beta);
[PnA,Pe]=vsc_func3_new(E,N1,N2,bndry,R,Pe,dVdt,para);
Pn=[PnA;Pnv];
[dVdtv,Fin,Fout]=vsc_func4(E,N1,N2,bndry,R,Pn,Pe);
[dSdt,OE1,PO2_new]= vsc_func5_Avg_new(S,E,N1,N2,V,Fin,Fout,L0,newd,PO2tisv,sen);
[dCMRO2ndt,OEn]=vsc_func7_Avg_bsk_mod_final(OE1,bndry,neu,CMRO2nref,W_gauss,file_xy,CMRO2n,PO2tis,PO2tismax);
[dATPdt,Jatp]=vsc_func9_Avg_new_final(CMRO2n,ATP,neu,lissom,CMRO2nref,ATPref);
dPO2tisdt=vsc_func8_Avg(OEn,CMRO2n,L0,bndry,nE);
dVdt=dVdtv;
dxdt=[dVdt;dSdt;dPO2tisdt;dATPdt;dbetadt;dCMRO2ndt];
if numel(find(isnan(dxdt))>0)>0
     disp('warning dxdt has gone wrong');
end
disp(t)

