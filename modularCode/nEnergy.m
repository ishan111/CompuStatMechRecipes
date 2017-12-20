if attpConf ~= 1
delpConf = (attpConf(pSelect,:) - attpConf);
delpConf = abs(delpConf);
delpConf_ = abs (delpConf-sysLen);
delpConf = min(delpConf, delpConf_) ;
delpConf = delpConf' ;
disp = dot(delpConf,delpConf).^0.5;
disp(pSelect)=[];
nE = sum(Interaction(disp));
else
    nE = 0 ;
end

