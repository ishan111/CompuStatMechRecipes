% HD simulation
mu = -1;
Kb= 1;
lambda = 1;
T = 1;
diaHD = 1;
len = [50 50 50];
dimes = length(len) ;
boundsNoP = [1 2500];
initNoP = 400 ;
noTrials = 1000000;
moveProbs=[1 2];
moveProbs = moveProbs./sum(moveProbs) ;
moveProbs= cumsum(moveProbs);
noMoves = length(moveProbs);
mExch = 2 ;
mDisp = 1 ;
mTot = noMoves+1;
dispLen = 20 ;
nAtt = zeros(1,mTot);
nAcc = zeros(1,mTot) ;
% accProb_l = mu/(Kb*T)- 3*log(lambda);
zz = .0115 ;
accProb_l = log(zz);
readConf = true;
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
pNoHist=zeros(1,boundsNoP(2)-boundsNoP(1)+1);
prodRun = true ;

tmmcC = zeros(3,boundsNoP(2)-boundsNoP(1)+1);
tmmcN = zeros(3,boundsNoP(2)-boundsNoP(1)+1);
tmmcUpdateFreq = 10000 ;
tmmcBias = false ;

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
            npConf(pSelect,:)= mod(npConf(pSelect,:),len);
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
                    accProb1 = min(exp(accProb),1);
                    tmmcC(1,pNo-boundsNoP(1)+1) =...
                        tmmcC(1,pNo-boundsNoP(1)+1) + accProb1 ;
                    tmmcC(2,pNo-boundsNoP(1)+1) =...
                        tmmcC(2,pNo-boundsNoP(1)+1) + 1 - accProb1 ;
                    if tmmcBias && trialNo > 500000 && tmmcN(1,pNo-boundsNoP(1)+1)~=0
                        bias = tmmcN(3,npNo-boundsNoP(1)+1)/tmmcN(1,pNo-boundsNoP(1)+1);
                        bias = log(bias);
                    else
                        bias = 0 ;
                    end
                    biasedAccProb = bias + accProb;
                    if log(rand())<biasedAccProb
                        nAcc(move) = nAcc(move)+1;
                        nAcc(mTot) = nAcc(mTot)+1 ;
                        pConf = npConf ;
                       
                       
                        pNo=npNo ;
                    end
                    if prodRun
                        pNoHist(pNo-boundsNoP(1)+1)= pNoHist(pNo-boundsNoP(1)+1)+1; 
                    end
                end
            elseif npNo > boundsNoP(1) %removal
                pSelect = randi(npNo);
                npConf(pSelect,:)=[] ;
                npNo = npNo - 1;
                %exp(beta*mu)/lambda^3
                %N/(pow(L,3) * zz)
                accProb = - accProb_l - log(prod(len)) + log(npNo +1) ;
                accProb1 = min(exp(accProb),1);
                tmmcC(3,pNo-boundsNoP(1)+1) =...
                    tmmcC(3,pNo-boundsNoP(1)+1) + accProb1 ;
                tmmcC(2,pNo-boundsNoP(1)+1) =...
                    tmmcC(2,pNo-boundsNoP(1)+1) + 1 - accProb1 ;
                % bias = P(t->tn)/P(tn->t)
                if tmmcBias && trialNo > 50000 && tmmcN(3,pNo-boundsNoP(1)+1)~=0
                    bias = tmmcN(1,npNo-boundsNoP(1)+1)/tmmcN(3,pNo-boundsNoP(1)+1);
                    bias = log(bias);
                else
                    bias = 0;
                end
                    biasedAccProb = bias + accProb ;
                
               if log(rand()) < biasedAccProb
                   nAcc(move) = nAcc(move)+1;
                   nAcc(mTot) = nAcc(mTot)+1 ;
                   pConf = npConf ;
                   
                   pNo=npNo ;
                   
               end
               if prodRun
                   pNoHist(pNo-boundsNoP(1)+1)= pNoHist(pNo-boundsNoP(1)+1)+1;
               end
            end
           
    end
    if mod(trialNo,tmmcUpdateFreq)
        tmmcN=tmmcC./sum(tmmcC);
    end
end
if dimes == 2
    scatter(pConf(:,1),pConf(:,2));
elseif dimes == 3
    scatter3(pConf(:,1),pConf(:,2),pConf(:,3));
end
figure(2)
 plot(boundsNoP(1):1:boundsNoP(2),pNoHist);
figure(3)
tmmcHist = ones(1,length(tmmcN(1,:)));
for i=2:length(tmmcN(1,:))
    tmmcHist(i)=tmmcN(1,i)/tmmcN(3,i-1);
end
tmmcHist(isnan(tmmcHist))= 1;
tmmcHist(isinf(tmmcHist))=1;
tmmcHist = cumprod(tmmcHist);
plot(log(tmmcHist))