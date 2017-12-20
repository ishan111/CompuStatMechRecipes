tmmcInit=true;
readConf = false;
tmmcBias = true ;

mu = 1;
Kb= 1; %1.38064852 × 10-23 m2 kg s-2 K-1
lambda = 1;
T = 1.2 ;
zz=0.055;
diaHD = 1;
sysLen = [8 8 8];
Vol=prod(sysLen);
dimes = length(sysLen) ;
pNbounds = [1 800];
initNoP = 100 ;
eps = 1 ;
sig = 1 ;
Interaction = @(dist) (4*eps.*((sig./dist).^12-(sig./dist).^6));

Ntrials = 100000;
moveProb=[1 3];
moveProb = moveProb./sum(moveProb) ;
moveProb= cumsum(moveProb);
noMoves = length(moveProb);
mExch = 2 ;
mDisp = 1 ;
mIns= 3 ;
mRem = 4 ;
mTot = noMoves+3;
maxDispLen = 3 ;
Natt = zeros(1,mTot);
Nacc = zeros(1,mTot) ;
% accProb_l = -E/(Kb*T) + mu/(Kb*T)- 3*log(lambda);
accProb_l = log(zz) ;

isCanon=false;
if prod(moveProb==[1 0])==1
    accProb_l=0;
    isCanon = true;
end

if ~readConf
    switch dimes
        case 1
            cpConf = linspace(0.1,sysLen-.1,initNoP) ;
        case 2
            cpConf = [linspace(0.1,sysLen(1)-.1,initNoP);linspace(0.1,sysLen(2)-.1,initNoP)]' ;
        case 3
            cpConf = [linspace(0.1,sysLen(1)-.1,initNoP);linspace(0.1,sysLen(2)-.1,initNoP);linspace(0.1,sysLen(3)-.1,initNoP)]';
    end
    cpNo = initNoP ;
end
pNoHist=zeros(1,pNbounds(2)-pNbounds(1)+1);
prodRun = true ;

if tmmcInit
    tmmcC = zeros(3,pNbounds(2)-pNbounds(1)+1);
    tmmcN = zeros(3,pNbounds(2)-pNbounds(1)+1);
end
% tmmcN=.055*ones(3,boundsNoP(2)-boundsNoP(1)+1);

tmmcUpdateFreq = 1000 ;
tmmcNupstart = 200000 ;
tmmcCupstart = 100000 ;
tmmcBiasStart= 210000;
a=zeros(1,Ntrials);
energy=zeros(1,Ntrials);
isProdRun = false;