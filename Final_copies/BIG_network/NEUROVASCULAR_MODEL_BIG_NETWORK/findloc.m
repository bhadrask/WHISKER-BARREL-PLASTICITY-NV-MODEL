function[loc_vas]=findloc(filename)
load(filename);
Gx=X(find(lev==4));
Gy=Y(find(lev==4));
G=[Gx Gy];
c=0;
caps=sqrt(length(Gx));
loc_vas=zeros(size(Gx));

for cj=1:(caps)
     for ci=1:(caps)
c=c+1;
cord=repmat([ci,cj],size(G,1),1);
diff_cord=G-cord;
loc_vas(c)=find(abs(diff_cord(:,1))+abs(diff_cord(:,2))==0);

     end
 end