m_tmn = m_tm ;
m_tmn(:,2:4) = m_tmn(:,2:4)./m_tm(:,2) ;
m_tm_hist = ones(1,length(m_tmn(1,:)));
for i=2:length(m_tmn(:,1))
    m_tm_hist(i)=m_tmn(i-1,3)/m_tmn(i,4);
end
 m_tm_hist(isnan(m_tm_hist))= 1;
 m_tm_hist(isinf(m_tm_hist))=1;
m_tm_hist = cumprod(m_tm_hist);
plot(m_tmn(:,1),log((m_tm_hist(:))./m_tm_hist(1)),'o-')
