e_tmn = e_tm ;
e_tmn(:,3:6) = e_tmn(:,3:6)./e_tm(:,2) ;
e_tm_hist = ones(1,length(e_tmn(1,:)));
for i=5:1:length(e_tmn(:,1))
    e_tm_hist(i)=e_tmn(i-1,3)/e_tmn(i,5);
end
% 
%  e_tm_hist(isnan(e_tm_hist))= 1;
%  e_tm_hist(isinf(e_tm_hist))=1;
e_tm_hist = cumprod(e_tm_hist);
plot(e_tmn(:,1),((e_tm_hist(:))./(e_tm_hist(7))))
