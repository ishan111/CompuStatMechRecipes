Natt(mTot) = Natt(mTot) + 1;
moveSelector = rand();
selectedMove = find(-moveSelector + moveProb>0,1);
attpConf = cpConf ;
attpNo = cpNo ;

switch(selectedMove)
    case(1)%disp
        Natt(mDisp)= Natt(mDisp) + 1;
        trial_disp;
    case(2)%exch
        Natt(mExch)= Natt(mExch) + 1;
        trial_exch;
end

