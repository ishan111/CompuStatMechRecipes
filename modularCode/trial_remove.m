Natt(mRem) = Natt(mRem) + 1 ;
if attpNo > pNbounds(1)
    pSelect=randi(attpNo);
    attpConf(pSelect,:)=[];
    attpNo = attpNo - 1;
    delPno = -1 ;
end