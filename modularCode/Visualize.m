figure(4)
cpNo
tmmcN=tmmcC./sum(tmmcC);
tmmcN(isnan(tmmcN)) =  0 ;
tmmcHist = zeros(1,length(tmmcN(1,:)));
for i=2:length(tmmcN(1,:))
if tmmcN(1,i-1)~=0 && tmmcN(3,i)~=0
tmmcHist(i)=log(tmmcN(1,i-1))-log(tmmcN(3,i));
else
tmmcHist(i) = 0;
end
end
% tmmcHist(isnan(tmmcHist))= 0;
% tmmcHist(isinf(tmmcHist))=0;
tmmcHist = cumsum(tmmcHist);
plot(tmmcHist,'g')
figure(1)
hold on;
plot(pNbounds(1):1:pNbounds(2),log(pNoHist),'r');