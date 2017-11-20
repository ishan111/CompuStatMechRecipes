figure(20);
clear a;
tmmc=log(transitionsProbMat);
tmmc(isinf(tmmc))=0;
a=zeros(1,length(length(tmmc(1,:))));
for i=2:length(tmmc(1,:))
    a(i)=-tmmc(i-1,i)+tmmc(i,i-1);
end
b=cumsum(a);
plot(b(30:100))
