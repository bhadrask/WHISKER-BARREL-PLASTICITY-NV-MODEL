function[lissom]=find_OSI(lissom,some,label,test)
if test==1 || some>2
xxnorm=lissom.layers{some}.xxnorm;
else
    xxnorm=lissom.layers{some}.Avg_neural_sheet_resp;
end
for kk=1:size(xxnorm,3)
    lr=2*label(kk)*ones(size(xxnorm,1),size(xxnorm,2)); % 10 x 10 full of 180 135 90 45
    mp_cos(:,:,kk)=cosd(lr).*xxnorm(:,:,kk);                % 10 x 10 x 1 mp=lr*xx
    mp_sin(:,:,kk)=sind(lr).*xxnorm(:,:,kk);
end
for pi=1:size(xxnorm,1)
    for pj=1:size(xxnorm,2)
        add_cos(pi,pj)=sum(mp_cos(pi,pj,:));
        add_sin(pi,pj)=sum(mp_sin(pi,pj,:));
        lissom.layers{some}.sigmap(pi,pj)=(atand(add_cos(pi,pj)./(add_sin(pi,pj))));  % orientation preference OP
        Sety(pi,pj)=sqrt(add_cos(pi,pj).^2 + add_sin(pi,pj).^2);
        add_xxnorm(pi,pj)=sum(xxnorm(pi,pj,:));
        if add_xxnorm(pi,pj)==0
            lissom.layers{some}.Selectivity(pi,pj)=0;
        else
            lissom.layers{some}.Selectivity(pi,pj)=Sety(pi,pj)./add_xxnorm(pi,pj);
        end
            % orientation selectivity
        if lissom.layers{some}.sigmap(pi,pj)<=0
            lissom.layers{some}.sigmap(pi,pj)=180 + lissom.layers{some}.sigmap(pi,pj);              % orientation preference
        end
    end
end