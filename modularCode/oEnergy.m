if cpConf ~= 1
    delpConf = (cpConf(pSelect,:) - cpConf);
    delpConf = abs(delpConf);
    delpConf_ = abs (delpConf-sysLen);
    delpConf = min(delpConf, delpConf_) ;
    delpConf = delpConf' ;
    disp = dot(delpConf,delpConf).^0.5;
    disp(pSelect)=[];
    oE = sum(Interaction(disp));
else
    oE = 0 ;
end
