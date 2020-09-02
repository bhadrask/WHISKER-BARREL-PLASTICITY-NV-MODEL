function[threshold]= fnctn(thr,ATP)
diff=ATP-thr;
threshold=zeros(size(ATP));
for i=1:size(ATP,1)
    for j=1:size(ATP,2)
       
if diff(i,j)>-0.1
    threshold(i,j)=0;
else
    threshold(i,j)=5;
end
    end
end