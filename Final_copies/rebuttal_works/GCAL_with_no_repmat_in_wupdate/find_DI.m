function[lissom]=find_DI(lissom,some,q,gg,test)

if test==1|| some>2
xxnorm=lissom.layers{some}.xxnorm;
else
    xxnorm=lissom.layers{some}.Avg_neural_sheet_resp;
end
for num_x=1:size(xxnorm,1)
    for num_y=1:size(xxnorm,2)
        if q(num_x,num_y)<=4
            q_opp(num_x,num_y)=q(num_x,num_y)+4;
            gg_opp(num_x,num_y)=xxnorm(num_x,num_y,q_opp(num_x,num_y));
            lissom.layers{some}.di(num_x,num_y)=1-gg_opp(num_x,num_y)/gg(num_x,num_y);
        end
        if q(num_x,num_y)<=8 && q(num_x,num_y)>4
            q_opp(num_x,num_y)=q(num_x,num_y)-4;
            gg_opp(num_x,num_y)=xxnorm(num_x,num_y,q_opp(num_x,num_y));
            lissom.layers{some}.di(num_x,num_y)=1-gg_opp(num_x,num_y)/gg(num_x,num_y);
        end
    end
end