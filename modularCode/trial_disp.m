pSelect = randi(cpNo) ;
dispLen = maxDispLen*rand() ;
dispVec = randn(1,dimes) ;
dispVec = dispVec/dot(dispVec,dispVec) ;
dispVec = dispLen*dispVec ;
attpConf(pSelect,:) = attpConf(pSelect,:) + dispVec ;
attpConf(pSelect,:) = mod(attpConf(pSelect,:), sysLen) ;

