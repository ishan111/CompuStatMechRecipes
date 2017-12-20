clear delpConf ;
if selectedMove == mRem
    nE=0;
else
    nEnergy;
end
if selectedMove == mIns
    oE=0;
else
    oEnergy;
end
if selectedMove == mRem && cpNo==pNbounds(1)
    oE = 0 ;
end
if selectedMove == mIns && cpNo==pNbounds(2)
    nE = 0 ;
end
    
dE = nE - oE ;
