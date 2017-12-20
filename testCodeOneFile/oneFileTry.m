T=1.2;    kb=1;
sysLen=[8 8 8];
zz = .055;
initpNo = 100 ;
noTrials = 1000000;
dimes=length(sysLen);
pNbounds = [2 200];
moveProb = [1 3] ;
maxDispLen = 3 ;
sig = 1;
eps = 1;
Interaction = @(r)((4.*eps).*(((sig./r).^12)-((sig./r).^6)));

readConf = false ;

%initialize
Vol = prod(sysLen);
lnvol = log(Vol);
accProb_pNo = log(zz);
pNo=initpNo;
pConf = zeros(pNo,3);
if ~readConf
    for i = 1:dimes
        pConf(:,i) = linspace(.1,sysLen(i)-.1,pNo);
    end
end

moveProb = moveProb/sum(moveProb);
moveProb = cumsum(moveProb);
mDisp = 1;
mExch = 2;
mTot = length(moveProb)+1 ;

pNoHist = zeros(1,pNbounds(2)-pNbounds(1)+1);

%counters
Natt = zeros(1,mTot);
Nacc = zeros(1,mTot);

energy = zeros(1,noTrials);


for trialNo=1:noTrials
    Natt(mTot) = Natt(mTot) + 1 ;
    npNo = pNo ;
    npConf = pConf ;
    mSelectRand = rand();
    mSelect = find(moveProb>mSelectRand,1);
    Natt(mSelect) = Natt(mSelect) + 1 ;
    switch mSelect
        case mDisp
            
            pSelect = randi(pNo);
            disp = randn(1,dimes);
            dispLen = maxDispLen*rand;
            disp = disp./dot(disp,disp);
            disp = disp*dispLen;
            npConf(pSelect,:) = npConf(pSelect,:) + disp ;
            npConf(pSelect,:) = mod(npConf(pSelect,:),sysLen);
            
            ndist_vector = npConf-npConf(pSelect,:) ;
            ndist_vector = abs(ndist_vector) ;
            ndist_vector_ = abs(ndist_vector - sysLen) ;
            ndist_vector = min(ndist_vector_,ndist_vector) ;
            ndist_vector = ndist_vector' ;
            ndist = dot(ndist_vector,ndist_vector).^0.5 ;
            ndist(pSelect)=[];
            nE = sum(Interaction(ndist));
            
            
            odist_vector = pConf-pConf(pSelect,:) ;
            odist_vector = abs(odist_vector) ;
            odist_vector_ = abs(odist_vector - sysLen) ;
            odist_vector = min(odist_vector_,odist_vector) ;
            odist_vector = odist_vector' ;
            odist = dot(odist_vector,odist_vector).^0.5 ;
            odist(pSelect)=[];
            oE = sum(Interaction(odist));
            
            dE = nE-oE ;
            accProb_e = -dE/(kb*T);
            accProb_L = accProb_e ;
            
            if log(rand)<=accProb_L
                pNo = npNo ;
                pConf = npConf ;
                Nacc(mSelect) = Nacc(mSelect) + 1 ;
                Nacc(mTot) = Nacc(mTot) + 1 ;
            end
            
        case mExch
            
            if rand() < 0.5 && pNo < pNbounds(2) %insertion
                npNo = npNo + 1;
                pSelect = npNo;
                npConf(pSelect,:) = rand(1,3).*sysLen ;
                
                ndist_vector = npConf-npConf(pSelect,:) ;
                ndist_vector = abs(ndist_vector) ;
                ndist_vector_ = abs(ndist_vector - sysLen) ;
                ndist_vector = min(ndist_vector_,ndist_vector) ;
                ndist_vector = ndist_vector' ;
                ndist = dot(ndist_vector,ndist_vector).^0.5 ;
                ndist(pSelect)=[];
                nE = sum(Interaction(ndist));
                
                
                oE = 0 ;
                
                dE = nE-oE ;
                accProb_e = -dE/(kb*T);
                accProb_L = accProb_e+accProb_pNo+lnvol-log(npNo) ;
                
                if log(rand)<=accProb_L
                    pNo = npNo ;
                    pConf = npConf ;
                    Nacc(mSelect) = Nacc(mSelect) + 1 ;
                    Nacc(mTot) = Nacc(mTot) + 1 ;
                    
                end
            elseif pNo > pNbounds(1)
                
                pSelect = randi(pNo);
                npConf(pSelect,:) = [] ;
                npNo = npNo - 1 ;
                
                nE = 0 ;
                
                
                odist_vector = pConf-pConf(pSelect,:) ;
                odist_vector = abs(odist_vector) ;
                odist_vector_ = abs(odist_vector - sysLen) ;
                odist_vector = min(odist_vector_,odist_vector) ;
                odist_vector = odist_vector' ;
                odist = dot(odist_vector,odist_vector).^0.5 ;
                odist(pSelect)=[];
                oE = sum(Interaction(odist));
                
                dE = nE-oE ;
                accProb_e = -dE/(kb*T);
                accProb_L = accProb_e-accProb_pNo-lnvol+log(npNo+1) ;
                
                if log(rand)<=accProb_L
                    pNo = npNo ;
                    pConf = npConf ;
                    Nacc(mSelect) = Nacc(mSelect) + 1 ;
                    Nacc(mTot) = Nacc(mTot) + 1 ;
                    
                end
            end
    end
       pNoHist(pNo-pNbounds(1)+1) = pNoHist(pNo-pNbounds(1)+1) + 1 ;
       energy(trialNo) = energy(trialNo) + dE ;
end
avEnergy = cumsum(energy)./(1:noTrials) ;

plot(log(pNoHist))
figure(2)
plot(avEnergy)

