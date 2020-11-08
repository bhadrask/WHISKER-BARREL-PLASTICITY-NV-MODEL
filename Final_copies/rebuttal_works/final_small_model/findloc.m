function[loc_vas]=findloc(filename,file_vessel_chara)
load(filename);
load(file_vessel_chara,'layers');
cap_layer=(numel(layers)+1)/2;
Gx=X(find(lev==cap_layer));
Gy=Y(find(lev==cap_layer));
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