delEnergy;
switch(selectedMove)
    case(mDisp)
        accProb_L = -dE/(Kb*T);
    case(mIns)
        accProb_L = -(dE/(Kb*T)) + log(Vol) - log (attpNo) + log(zz)  ;
    case(mRem)
        accProb_L = -(dE/(Kb*T)) - log(Vol) + log(cpNo) - log(zz);
end