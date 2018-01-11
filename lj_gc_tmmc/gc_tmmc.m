T=1.2;
% Mu = -3.9 ;
% zz = exp(Mu/T) ;
zz=0.052;
sysLen=[8 8 8] ;
Nmin = 2;
Nmax = 500 ;

Ntrials = 1E7 ;
equilRunNo=1E5;

pNo=20;

moveProb = [.25 1] ;
maxDisp = 3 ;

dimes = length(sysLen);
V = prod(sysLen);

sig = 1;
eps = 1 ;
Interaction = @(dist) (4*eps.*(((sig./dist).^12)-((sig./dist).^6)));
pNoHist = zeros(1,Nmax-Nmin+1);
tmmcC = zeros(Nmax-Nmin+1,3) ;

pConf= zeros(pNo,3);
for i=1:dimes
    pConf(:,i)=linspace(0.1,sysLen(i)-0.1,pNo);
end


for trialNo = 1:Ntrials
    move = select_move(moveProb);
    switch(move)
        case(1) %displace
            pSelect = randi(pNo);
            disp_vec = randn(1,dimes);
            disp_vec = disp_vec./sum(disp_vec.^2)^0.5 ;
            disp = maxDisp*disp_vec ;
            newPos = pConf(pSelect,:) + disp ;
            newPos = mod(newPos,sysLen);
            dist = dist_calc(dimes,sysLen,pSelect,pConf,newPos) ;
            dE = sum(Interaction(dist));
            dist = dist_calc(dimes,sysLen,pSelect,pConf,pConf(pSelect,:)) ;
            dE = dE - sum(Interaction(dist)) ;
            acc_prob_l = acc_prob_ln_disp(T,dE) ;
            if log(rand())<acc_prob_l
                pConf(pSelect,:)=newPos;
            end
            
            
        case(2) %exchange
            
            if rand()<0.5 %insert
                if pNo<Nmax
                    
                    newPos = sysLen.*rand(1,3) ;
                    dist = dist_calc(dimes,sysLen,0,pConf,newPos) ;
                    dE = sum(Interaction(dist));
                    acc_prob_l = acc_prob_ln_ins( T,zz,V,pNo,dE ) ;
                    if log(rand())<acc_prob_l
                        pNo = pNo+1;
                        pConf(pNo,:)=newPos;
                    end
                    acc_prob = exp(acc_prob_l);
                    if trialNo>equilRunNo
                        tmmcC(pNo-Nmin+1,3) = tmmcC(pNo-Nmin+1,3) + acc_prob;
                        tmmcC(pNo-Nmin+1,2) = tmmcC(pNo-Nmin+1,2) + 1 - acc_prob;
                    end
                else
                    newPos = sysLen.*rand(1,3) ;
                    dist = dist_calc(dimes,sysLen,0,pConf,newPos) ;
                    dE = sum(Interaction(dist));
                    acc_prob_l = acc_prob_ln_ins( T,zz,V,pNo,dE ) ;
                    
                    acc_prob = exp(acc_prob_l);
                    if trialNo>equilRunNo
                        tmmcC(pNo-Nmin+1,3) = tmmcC(pNo-Nmin+1,3) + acc_prob;
                        tmmcC(pNo-Nmin+1,2) = tmmcC(pNo-Nmin+1,2) + 1 - acc_prob;
                    end
                end
            else   %remove
                if pNo>Nmin
                    pSelect = randi(pNo);
                    dist = dist_calc(dimes,sysLen,pSelect,pConf,pConf(pSelect,:)) ;
                    dE = - sum(Interaction(dist)) ;
                    acc_prob_l = acc_prob_ln_rem( T,zz,V,pNo,dE ) ;
                    if log(rand())<acc_prob_l
                        pNo=pNo-1;
                        pConf(pSelect,:)= [];
                        
                    end
                    acc_prob = exp(acc_prob_l);
                    if trialNo>equilRunNo
                        tmmcC(pNo-Nmin+1,1) = tmmcC(pNo-Nmin+1,1) + acc_prob;
                        tmmcC(pNo-Nmin+1,2) = tmmcC(pNo-Nmin+1,2) + 1 - acc_prob;
                    end
                    
                else
                    pSelect = randi(pNo);
                    dist = dist_calc(dimes,sysLen,pSelect,pConf,pConf(pSelect,:)) ;
                    dE = - sum(Interaction(dist)) ;
                    acc_prob_l = acc_prob_ln_rem( T,zz,V,pNo,dE ) ;
                    acc_prob = exp(acc_prob_l);
                    if trialNo>equilRunNo
                        tmmcC(pNo-Nmin+1,1) = tmmcC(pNo-Nmin+1,1) + acc_prob;
                        tmmcC(pNo-Nmin+1,2) = tmmcC(pNo-Nmin+1,2) + 1 - acc_prob;
                    end
                end
            end
            
    end
    if trialNo>equilRunNo
        pNoHist(pNo-Nmin+1) = pNoHist(pNo-Nmin+1)+1 ;
    end
end
hold on
for i=2:1-Nmin+Nmax
    tmmcHist(i)=-log(tmmcN(i-1,3))+log(tmmcN(i,1));
end
tmmcHist(isnan(tmmcHist))=0;
tmmcHist(isinf(tmmcHist))=0;
tmmcHist=cumsum(tmmcHist);
plot(tmmcHist-tmmcHist(find(tmmcHist,1)))
pNoHist_l=log(pNoHist);
pNoHist_l(isnan(pNoHist_l))=0;
pNoHist_l(isinf(pNoHist_l))=0;
plot(pNoHist_l-pNoHist_l(find(pNoHist_l,1)))
