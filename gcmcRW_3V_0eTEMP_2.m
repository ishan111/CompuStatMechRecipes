Nmin = 1 ;
Nmax = 700 ;
% Neighbor List was implemented
lambda = 1;
eps = 1 ;
sig = 1 ;
Interaction = @(r)(4*eps*((sig./r).^12-(sig./r).^6));
%% SYSTEM PARAMS
Kb = 1;
%1.38064852e-23 ;
T = .7 ;
Mu= -3;
Nparticles = 200 ;
Prob_Fact_Log = @ (E,N)(( Mu*N - E )/( Kb * T )); 
SysLen = [22 22];
Dimes = length(SysLen);
%% MOVE PARAMS
Ntrials = 100000000;
moves_vector = [1 2 2] ;% to be modified later 
nMoveType = length(moves_vector) ;
moves_prob_vector = moves_vector/sum(moves_vector) ;
moves_prob_vector = cumsum(moves_prob_vector) ; %Tower selector
maxdispLen = 15 ;
Nproperties = 2 ;
%% INTERACTION CALC CHOICES
isNlist = false ;
isClist = false ;
isLongRange = false ;
%% INTERACTION CALC PARAMS
Rc = 2.5 * sig ;
Rv = 2.7 * sig ;


%% INITIAL STATE
switch(Dimes)
    case(1)
        init_positions = linspace(0,SysLen(1)-.01,Nparticles);
    case(2)
        init_positions = [linspace(0,SysLen(1)-.01,Nparticles);...
           linspace(0,SysLen(2)-.01,Nparticles)];
    case(3)
        init_positions = [linspace(0,SysLen(1)-.01-.01,Nparticles);...
            linspace(0,SysLen(2)-.01,Nparticles);...
            linspace(0,SysLen(2)-.01,Nparticles)];

end
init_positions = init_positions';

trial(1,Nparticles) = struct('is_accepted',false,...
    'Nparticles',Nparticles,...
    'properties',zeros(1,Nproperties),...
    'SysLen',SysLen) ;

oconfig = struct('positions',init_positions,...
    'Nparticles',Nparticles,...
    'properties',zeros(1,Nproperties),...
    'SysLen',SysLen) ;

% equis_config_maker = linspace(0,SYSTEM_SIZE,nthroot(NO_PARTICLES,3)) ;
% [X_init_positions, Y_init_positions, Z_init_positions] = meshgrid(equis_config_maker,equis_config_maker,equis_config_maker) ;
%% MOVE INDICES
mTotal = nMoveType + 1 ;
mDisplace = 1 ;
mRemove = 2 ;
mInsert = 3 ;
%% ACCEPTANCE LOGS
Natt = zeros(1,nMoveType + 1) ;
Nacc = zeros(1,nMoveType + 1) ;
%% PROPERTY INDICES
pAv_energy = 2;
pEnergy = 1;
pDSn = zeros(1,700);
pavDSn = zeros(1,Ntrials) ;
ProdRunNo = 10000;
if isClist
    nMesh = 8 * ones(1,Dimes) ; %#ok<*UNRCH>
end
%% NEIGHBOR LIST
if isNlist
    Nlist(1,Nparticles) =  struct('neighbors',[],'nbrIndices',[]);

    for i = 1 : oconfig.Nparticles
        particleIndices = 1:1:oconfig.Nparticles ;
        dist_norm = DistCalc( oconfig.positions, i,oconfig.Nparticles, oconfig ) ;
        Nlist(i).neighbors = oconfig.positions(dist_norm<=Rv,:);
        Nlist(i).nbrIndices = particleIndices(dist_norm<=Rv) ;
    end
end
%% CELL LIST
% if isClist
%     Clist(nMesh) = struct('nParticles',0,'positions',[]);
%     switch Dimes
%         case 1
%             gMesh = linspace(0, SysLen, nMesh + 1) ;
%             for meshX=2:nmesh+1
%                 Clist(meshX-1).nParticles = 0 ;
%             end
%              
%             for particleNo = 1:oconfig.Nparticles
%                 for meshX=2:nmesh+1
%                     if oconfig.positions(particleNo)<=gMesh(meshx)
%                        
%                         Clist(meshX-1).nParticles = Clist(meshX-1).nParticles + 1;
%                         Clist(meshX-1).positions(Clist(meshX-1).nParticles)...
%                             = oconfig.positions(particleNo);
%                         break;
%                         
%                     end
%                 end
%             end
%         case 2
%             gMesh = [linspace(0, SysLen(1), nMesh(1) + 1);...
%                 linspace(0, SysLen(2), nMesh(2) + 1)] ;
%             for meshX=2:nMesh(1)+1
%                 for meshY = 2:nMesh(2)+1
%                     
%                     Clist(meshX-1,meshY-1).nParticles = 0 ;
%                     
%                 end
%             end
%             
%             for particleNo = 1:oconfig.Nparticles
%                 bFlag = false ;
%                 for meshX=2:nMesh(1)+1
%                     for meshY = 2:nMesh(2)+1
%                         if (oconfig.positions(particleNo,1)<=gMesh(1,meshX))
%                             if(oconfig.positions(particleNo,2)<=gMesh(2,meshY))
%                                 Clist(meshX-1,meshY-1).nParticles = Clist(meshX-1,meshY-1).nParticles + 1;
%                                 %                             Clist(meshX-1,meshY-1).positions(Clist(meshX-1,meshY-1).nParticles,:)...
%                                 %                                 = oconfig.positions(particleNo,:);
%                                 bFlag = true ;
%                                 break;
%                                 
%                             end
%                             
%                         end
%                        
%                     end
%                     if bFlag==true
%                         break
%                     end
%                 end
%             end
%         case 3
%             gMesh = [linspace(0, SysLen(1), nMesh(1) + 1);...
%                 linspace(0, SysLen(2), nMesh(2) + 1);...
%                 linspace(0, SysLen(3), nMesh(3) + 1)] ;
%             for meshX=2:nMesh(1)+1
%                 for meshY = 2:nMesh(2)+1
%                     for meshZ =  2:nMesh(2)+1
%                         
%                         Clist(meshX-1,meshY-1,meshZ-1).nParticles = 0 ;
%                     end
%                 end
%             end
%            
%             for particleNo = 1:oconfig.Nparticles
%                  bFlag = false ;
%                 for meshX=2:nMesh(1)+1
%                     for meshY = 2:nMesh(2)+1
%                         for meshZ =  2:nMesh(2)+1
%                             if oconfig.positions(particleNo,1)<=gMesh(1,meshX)&&...
%                                    oconfig.positions(particleNo,2)<=gMesh(2,meshX)&&...
%                                   oconfig.positions(particleNo,3)<=gMesh(3,meshZ)
%                                 Clist(meshX-1,meshY-1,meshZ-1).nParticles = Clist(meshX-1,meshY-1,meshZ-1).nParticles + 1;
%                                 Clist(meshX-1,meshY-1,meshZ-1).positions(Clist(meshX-1,meshY-1,meshZ-1).nParticles)...
%                                     = oconfig.positions(particleNo);
%                                 bFlag = true ;
%                                 break;
%                                 
%                             end
%                             if bFlag == true
%                                 break;
%                             end
%                             if bFlag == true
%                                 break
%                             end
%                         end
%                     end
%                 end
%             end
%             
%     end
%     
% end
%% TMMC PARAMS
isTMMC = true ;
isTMMCbias = true ;
isTMMC2 = false ;
TMmeshNo = 10*ones(size(SysLen));
if isTMMC2
    rhoPos = [4 , 5 ; 2 , 3];
    for i=1:Dimes
        tmmcMesh(:,i)= linspace(0,SysLen(i),TMmeshNo(i)+1); 
    end
    tmmcRho = zeros(TMmeshNo) ;
    

    switch Dimes
        case 1
            for i= 1:oconfig.Nparticles
                meshP = floor(oconfig.positions(i,:)./(SysLen./TMmeshNo)) + 1;
                
                tmmcRho(meshP) = tmmcRho(meshP) + 1 ;
            end
        case 2
            for i= 1:oconfig.Nparticles
                meshP = floor(oconfig.positions(i,:)./(SysLen./TMmeshNo)) + 1;
                
                tmmcRho(meshP(1),meshP(2)) = tmmcRho(meshP(1),meshP(2)) + 1 ;
            end
        case 3
            for i= 1:oconfig.Nparticles
                meshP = floor(oconfig.positions(i,:)./(SysLen./TMmeshNo)) + 1;
                
                tmmcRho(meshP(1),meshP(2),mesh(3)) = tmmcRho(meshP(1),meshP(2),mesh(3)) + 1 ;
            end
            
    end
    transitionsMatrix = zeros(11,11,11,11) ;
    transitionsProbMat = transitionsMatrix ;
end
if isTMMC == true
    
    %    tmMacroProp ;
    transitionsMatrix = zeros(700,700) ;
    transitionsProbMat = transitionsMatrix ;
end

%% THE MAIN LOOP (SAMPLE GEN)
for trialNo = 1 : Ntrials
    Natt(mTotal) = Natt(mTotal) + 1 ;
    trial_config = oconfig ;
    if isNlist
        trial_Nlist = Nlist ;
    end
    %% MAKING A MOVE
    movselector = rand();
    for i = 1:nMoveType
        if movselector <= moves_prob_vector(i)
            moveSelect = i ;
            break;
        end
    end
    %     moveSelect = randi(nMoveType);
    Natt(moveSelect ) = Natt(moveSelect) + 1 ;
    
    switch(moveSelect)
        case(mDisplace)
            %             disp = rand(1,Dimes) ;
            disp = normrnd(0,1,[1 Dimes+2]);
            %             dispLen = rand() * maxdispLen ;
            dispLen = maxdispLen;
            disp = (dispLen/dot(disp,disp)) * disp ;
            pSelect = randi(trial_config.Nparticles) ;
            new_pos = trial_config.positions(pSelect) ;
            new_pos= new_pos + disp(1:Dimes) ;
            new_pos = mod( new_pos, trial_config.SysLen );
            trial_config.positions(pSelect,:) = new_pos ;
            if isNlist
                if dispLen <= Rv
                    for i = 1 : trial_config.Nparticles
                        particleIndices = 1:1:trial_config.Nparticles ;
                        dist_norm = DistCalc( trial_config.positions, i,trial_config.Nparticles, trial_config ) ;
                        trial_Nlist(i).neighbors = trial_config.positions(dist_norm<=Rv,:);
                        trial_Nlist(i).nbrIndices = particleIndices(dist_norm<=Rv) ;
                    end
                end
                
                %% OLD ENERGY CALC
                
                ocalcPos = Nlist(pSelect).neighbors;
                ocalcNparticles = length(ocalcPos(:,1));
                ocalcNparticles = ocalcNparticles + 1 ;
                ocalcPos(ocalcNparticles ,:)= oconfig.positions(pSelect,:) ;
                opSelectCalc = ocalcNparticles ;
            else
                ocalcPos = oconfig.positions ;
                opSelectCalc = pSelect ;
                ocalcNparticles = oconfig.Nparticles ;
            end
            %% DISTANCES CALCULATION
            dist_norm = DistCalc( ocalcPos,opSelectCalc,ocalcNparticles,oconfig );
            
            %% NET ENERGY OF YET UNMOVED PARTICLE
            oconfig_energy = sum(Interaction(dist_norm)) ;
            %% NEW ENERGY CALC
            if ~isNlist
                calcPos = trial_config.positions ;
                pSelectCalc = pSelect ;
                calcNparticles = trial_config.Nparticles ;
            else
                calcPos = trial_Nlist(pSelect).neighbors;
                calcNparticles = length(calcPos(:,1));
                calcNparticles = calcNparticles + 1 ;
                calcPos(calcNparticles ,:)= trial_config.positions(pSelect,:) ;
                pSelectCalc = calcNparticles ;
            end
            %% DISTANCES CALCULATION
            dist_norm = DistCalc( calcPos,pSelectCalc,calcNparticles,trial_config ) ;
            %% NET ENERGY OF MOVED PARTICLE
            trial_config_energy = sum(Interaction(dist_norm)) ;
        case(mRemove)
            if trial_config.Nparticles ~= Nmin
                pSelect = randi(trial_config.Nparticles) ;
                trial_config.Nparticles = trial_config.Nparticles - 1 ;
                new_pos = [] ;
                trial_config.positions(pSelect,:) = [] ;
                if isNlist
                    for i = 1 : trial_config.Nparticles
                        particleIndices = 1:1:trial_config.Nparticles ;
                        dist_norm = DistCalc( trial_config.positions, i,trial_config.Nparticles, trial_config ) ;
                        trial_Nlist(i).neighbors = trial_config.positions(dist_norm<=Rv,:);
                        trial_Nlist(i).nbrIndices = particleIndices(dist_norm<=Rv) ;
                    end
                end
                %% OLD ENERGY CALC
                if ~isNlist
                    ocalcPos = oconfig.positions ;
                    opSelectCalc = pSelect ;
                    ocalcNparticles = oconfig.Nparticles ;
                else
                    
                    ocalcPos = Nlist(pSelect).neighbors;
                    ocalcNparticles = length(ocalcPos(:,1));
                    ocalcNparticles = ocalcNparticles + 1 ;
                    ocalcPos(ocalcNparticles ,:)= oconfig.positions(pSelect,:) ;
                    opSelectCalc = ocalcNparticles ;
                end
                %% DISTANCES CALCULATION
                dist_norm = DistCalc( ocalcPos,opSelectCalc,ocalcNparticles,oconfig );
                %% NET ENERGY OF REMOVED PARTICLE
                
                oconfig_energy = sum(Interaction(dist_norm)) ;
                %% NET ENERGY AFTER
                
                trial_config_energy = 0 ;
            else
                trial_config.positions = oconfig.positions ;
                trial_config_energy = 0 ;
                oconfig_energy = 0 ;
            end
        case(mInsert)
            if trial_config.Nparticles ~= Nmax
                new_pos = SysLen.*rand(1,Dimes);
                trial_config.Nparticles = trial_config.Nparticles + 1 ;
                pSelect = trial_config.Nparticles ;
                trial_config.positions(pSelect,:) = new_pos ;
                if isNlist
                    for i = 1 : trial_config.Nparticles
                        particleIndices = 1:1:trial_config.Nparticles ;
                        dist_norm = DistCalc( trial_config.positions, i,trial_config.Nparticles, trial_config ) ;
                        trial_Nlist(i).neighbors = trial_config.positions(dist_norm<=Rv,:);
                        trial_Nlist(i).nbrIndices = particleIndices(dist_norm<=Rv) ;
                    end
                end
                %% NET ENERGY BEFORE INSERT
                oconfig_energy = 0 ;
                %% NEW ENERGY CALC
                if ~isNlist
                    calcPos = trial_config.positions ;
                    pSelectCalc = pSelect ;
                    calcNparticles = trial_config.Nparticles ;
                else
                    calcPos = trial_Nlist(pSelect).neighbors;
                    calcNparticles = length(calcPos(:,1));
                    calcNparticles = calcNparticles + 1 ;
                    calcPos(calcNparticles ,:)= trial_config.positions(pSelect,:) ;
                    pSelectCalc = calcNparticles ;
                end
                %% DISTANCES CALCULATION
                dist_norm = DistCalc( calcPos,pSelectCalc,calcNparticles,trial_config ) ;
                %% NET ENERGY OF INSERTED PARTICLE
                trial_config_energy = sum(Interaction(dist_norm)) ;
            else
                trial_config.positions = oconfig.positions ;
                oconfig_energy = 0 ;
                trial_config_energy = 0 ;
            end
    end
  
    delta_energy = trial_config_energy - oconfig_energy ;
    
    accProbLog = 0 ;
    %% accProbLog
    switch moveSelect
        case 1
            accProbLog = -(1/Kb*T).*delta_energy ;
        case 2
            accProbLog = -(1/Kb*T).*(delta_energy + Mu) + 3*log(lambda) + log(oconfig.Nparticles) - log(prod(SysLen));
        case 3
            accProbLog = -(1/Kb*T).*(delta_energy - Mu) - 3*log(lambda) - log(trial_config.Nparticles) + log(prod(SysLen));
    end
    accProbLog = min(0,accProbLog) ;
    accRand = log(rand());
   
      
    %% TMMC
    % find the macrostate
    % make tm based on certain meshed(maybe naturally) microstate
    % update it regularly
    bias = 0 ;
    if isTMMC2
        for i=1:Dimes
            tmmcMesh(:,i)= linspace(0,SysLen(i),TMmeshNo(i)+1);
        end
        tmmcRho = zeros(TMmeshNo) ;
        switch Dimes
            case 1
                for i= 1:oconfig.Nparticles
                    meshP = floor(oconfig.positions(i,:)./(SysLen./TMmeshNo)) + 1;
                    tmmcRho(meshP) = tmmcRho(meshP) + 1 ;
                end
                tmMacroProp = [ tmmcRho(rhoPos(1)) tmmcRho(rhoPos(2))] ;
            case 2
                for i= 1:oconfig.Nparticles
                    meshP = floor(oconfig.positions(i,:)./(SysLen./TMmeshNo)) + 1;
                    tmmcRho(meshP(1),meshP(2)) = tmmcRho(meshP(1),meshP(2)) + 1 ;
                end
                
                tmMacroProp = [ tmmcRho(rhoPos(1,1),rhoPos(1,2)) tmmcRho(rhoPos(2,1),rhoPos(2,2)) ] ;
            case 3
                for i= 1:oconfig.Nparticles
                    meshP = floor(oconfig.positions(i,:)./(SysLen./TMmeshNo)) + 1;
                    tmmcRho(meshP(1),meshP(2),mesh(3)) = tmmcRho(meshP(1),meshP(2),mesh(3)) + 1 ;
                end
                
                tmMacroProp = [ tmmcRho(rhoPOs(1,1),rhoPos(1,2),rhoPos(1,3))  tmmcRho(rhoPOs(2,1),rhoPos(2,2),rhoPos(2,3))] ;
        end
    end
    if  isTMMC && trialNo> 10000 && (moveSelect == 3 || moveSelect == 2)
        tmMacroProp = [ oconfig.Nparticles  trial_config.Nparticles] ;
        % WHAT IF THE PROPERTY IS NOT CHANGED
        transitionsMatrix(tmMacroProp(1),tmMacroProp(2)) = ...
            transitionsMatrix(tmMacroProp(1),tmMacroProp(2))...
            + exp(accProbLog) ;
         transitionsMatrix(tmMacroProp(1),tmMacroProp(1)) = ...
            transitionsMatrix(tmMacroProp(1),tmMacroProp(1))...
            + 1 - exp(accProbLog) ;
        tmmcNormFactor = sum(transitionsMatrix(tmMacroProp(1),:)) ;
        if tmmcNormFactor ~= 0 && mod(trialNo,1000)==0
        transitionsProbMat(tmMacroProp(1),:) = transitionsMatrix(tmMacroProp(1),:)...
            ./ tmmcNormFactor ;
        end
        
        
    end
    if isTMMC && trialNo >200000 && isTMMCbias && (moveSelect == 3 || moveSelect == 2)
        bias = log(transitionsProbMat(tmMacroProp(2),tmMacroProp(1)))-...
            log(transitionsProbMat(tmMacroProp(1),tmMacroProp(2))) ;
%         if isnan(bias) || isinf(bias)
%             bias =  0 ;
%         end
        
    end
    if isTMMC2 == true
        
    end
     biasedAccProb = accProbLog + bias ; 
    %%
    if accRand<=biasedAccProb
        Nacc(mTotal) = Nacc(mTotal) + 1;
        Nacc(moveSelect) = Nacc(moveSelect) + 1;
        oconfig = trial_config ;
        if isNlist
            Nlist = trial_Nlist ;
        end
    end
    if trialNo>=ProdRunNo
        pDSn(oconfig.Nparticles) = pDSn(oconfig.Nparticles) + 1 ;
        pavDSn(trialNo)=pavDSn(trialNo-1);
        pavDSn(trialNo)= ((trialNo-ProdRunNo) * pavDSn(trialNo) + oconfig.Nparticles)/(trialNo-ProdRunNo+1)  ;
    end
end
Untitled7;
