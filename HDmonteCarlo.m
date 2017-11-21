% HD simulation
mu = -1;
Kb= 1;
lambda = 1;
T = 1;
diaHD = 1;
len = [20 20];
dimes = length(len) ;
boundsNoP = [10 50];
initNoP = 15 ;
noTrials = 1000000;
moveProbs=[1 2];
moveProbs = moveProbs./sum(moveProbs) ;
moveProbs= cumsum(moveProbs);
noMoves = length(moveProbs);
mExch = 2 ;
mDisp = 1 ;
mTot = noMoves+1;
dispLen = 10 ;
nAtt = zeros(1,mTot);
nAcc = zeros(1,mTot) ;
readConf = true;
accProb_l = mu/(Kb*T)- 3*log(lambda);
if ~readConf
    switch dimes
        case 1
            pConf = linspace(0,len,initNoP) ;
        case 2
            pConf = [linspace(0,len(1),initNoP);linspace(0,len(2),initNoP)]' ;
        case 3
            pConf = [linspace(0,len(1),initNoP);linspace(0,len(2),initNoP);linspace(0,len(3),initNoP)]';
    end
    pNo = initNoP ;
end
for trialNo = 1:noTrials
    mSelRand= rand();
    for moveNo = 1 : noMoves
        if mSelRand<=moveProbs(moveNo)
            move = moveNo ;
            break;
        end
         
    end
    nAtt(move)=nAtt(move)+1;
    nAtt(mTot)=nAtt(mTot)+1;
    npConf = pConf ;
    npNo = pNo;
    switch move
        case mDisp 
            
            disp = normrnd(0,1,[1 dimes+2]);
            disp = (1/(dot(disp,disp))^.5) * disp ;
            pSelect = randi(npNo) ;
            npConf(pSelect,:) = pConf(pSelect,:) + dispLen*disp(1:dimes);
            dist = DistCalc(npConf,pSelect,npNo,len);
            %dE = sum(Interaction(dist));
            if ~any(dist<diaHD)
                nAcc(move) = nAcc(move)+1;
                nAcc(mTot) = nAcc(mTot)+1 ;
                pConf = npConf ;
                pNo=npNo ;
                
            end
        case mExch
            if rand()<.5 && npNo < boundsNoP(2) %insertion
                npNo = npNo + 1;
                pSelect = npNo ;
                npConf(npNo,:) = rand(1,dimes).*len ; 
                dist = DistCalc(npConf,pSelect,npNo,len);
                if ~any(dist<diaHD)
                    %exp(beta*mu)/lambda^3
                    %(pow(L,3) * zz)/(N+1)
                    accProb = accProb_l + log(prod(len)) - log (npNo);
                    if log(rand())<accProb
                        nAcc(move) = nAcc(move)+1;
                        nAcc(mTot) = nAcc(mTot)+1 ;
                        pConf = npConf ;
                        pNo=npNo ;
                    end
                    
                end
            elseif npNo > boundsNoP(1) %removal
                pSelect = randi(npNo);
                npConf(pSelect,:)=[] ;
                npNo = npNo - 1;
                %exp(beta*mu)/lambda^3
                %N/(pow(L,3) * zz)
                accProb = - accProb_l - log(prod(len)) + npNo +1 ; 
               if log(rand())<accProb
                    nAcc(move) = nAcc(move)+1;
                    nAcc(mTot) = nAcc(mTot)+1 ;
                    pConf = npConf ;
                    pNo=npNo ;
                end
            end
    end
end
if dimes == 2
    scatter(pConf(:,1),pConf(:,2));
elseif dimes == 3
    scatter3(pConf(:,1),pConf(:,2),pConf(:,3));
end
pNo