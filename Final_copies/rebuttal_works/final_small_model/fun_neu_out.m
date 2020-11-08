function dxdt=fun_neu_out(t,x,Ae,lissom,opts,Amax,IN,Z_constraint,W_gauss,loc_vas,file_xy,file_vess_char)
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
ATPref=2.5;
ATP=ATPref*ones(lissom.layers{2}.dim(1)* lissom.layers{2}.dim(2),1);
if (t<=tstart) || (t>tstop)
    input=zeros(size(input1));
end

t1=t-tau;
if (t1<=tstart) || (t1>tstop)
    input2=zeros(size(input1));
else
    input2=input1;
end
% if t1>tstop
%     keyboard;
% end
% if (t>tstart) ||(t1>tstart)||(t>tstop)||(t1>tstop)
%     keyboard;
% end


%% Loading previous iteration values to find the increments
beta=x;

lissom1=lissom;
lissomref=lissom;
Zold_ref=lissomref.layers{2}.Zold;
[lissom.layers{2}.Zold,lissom1.layers{2}.Zold]=getGlobalx;
 [cnt,neuout,ntime]=getGlobalneu;
 %% if stepdown happens in ode
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
%      kcnt=floor(10*abs((abs(ntime(cnt)-t)-0.1)));
% 
%      for i=0:kcnt
     neu=inst_responsev1(lissom,opts,input,Z_constraint,ATP,Amax);
%      lissom.layers{2}.Zold=neu;
      cnt=cnt+1;
    neuout(:,:,cnt)=neu;
    ntime(cnt)=t;
     setGlobalneu(cnt,neuout,ntime);
%      end
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
% figure(1);subplot(211); imagesc(neu);title(num2str(t))
% subplot(212);imagesc(neu_prev);title(num2str(t));pause(0.01);
% % if t1<=tstart
% %     lissom1.layers{2}.Zold=Zold_ref;
% % end
% % % neu=inst_response(lissom,opts,input,Z_constraint,ATP,Amax);
% % % neu_prev=inst_response(lissom1,opts,input2,Z_constraint,ATP,Amax);
% if t1>tstart
  setGlobalx(neu,neu_prev);
% else
%     setGlobalx(neu,lissom1.layers{2}.Zold);
% end

% % if ismember(round(t,3),round(0:0.1:4,3))
%    
%    
%     cnt=cnt+1;
%     neuout(:,:,cnt)=neu;
%     ntime(cnt)=t;
%      setGlobalneu(cnt,neuout,ntime);
% % end
dbetadt=beta_diff(neu_prev, bndry,Ae,beta,W_gauss,file_xy,loc_vas,N1,N2);

dxdt=dbetadt;
if numel(find(isnan(dxdt))>0)>0
     disp('warning dxdt has gone wrong');
end
disp(t)

