if selectedMove == mDisp && cpNo == 1
    cpConf = attpConf;
elseif cpNo == pNbounds(1) && selectedMove == mRem
    0;
elseif cpNo == pNbounds(2) && selectedMove == mIns
    0;
else
    if log(rand()) < accProb_L
        accepted=true;
        Nacc(mTot) = Nacc(mTot) + 1;
        Nacc(selectedMove) = Nacc(selectedMove) + 1 ;
        if selectedMove == mRem || selectedMove == mIns
            Nacc(mExch)= Nacc(mExch) + 1 ;
        end
        cpConf=attpConf;
        cpNo=attpNo;
    end
end