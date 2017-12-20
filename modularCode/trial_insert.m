Natt(mIns) = Natt(mIns) + 1 ;
if attpNo<pNbounds(2)
    attpNo = attpNo + 1 ;
    pSelect = attpNo ;
    attpConf(pSelect,:) = rand(1,3).*sysLen;
    delpno = +1 ;
end