if selectedMove == mIns
    accProb = min(1,exp(accProb_L));
    tmmcC(2,cpNo)= tmmcC(2,cpNo) + 1 - accProb ;
    tmmcC(1,cpNo)= tmmcC(1,cpNo) + accProb ;
end
    
if selectedMove == mRem
    accProb = min(1,exp(accProb_L));
    tmmcC(2,cpNo)= tmmcC(2,cpNo) + 1 - accProb ;
    tmmcC(3,cpNo)= tmmcC(3,cpNo) + accProb ;
end