%made by Ishan Arora under the guidance of  Dr.Rane
%We use 6-12 potentials,
%optimization tactics not employed yet
%also need to improve readability of the code
%replace by suitable(reduced)units later.

%No of Attempted and accepted moves

attempted_moves = 10000;
attempted_remove = 0 ;
attempted_add = 0 ;
attempted_disp = 0 ;
accepted_moves = 0;
accepted_add = 0 ;
accepted_remove = 0 ;
accepted_disp = 0 ;

%Markov Chain Moves characteristics

moves_vector = [1 1 1] ;
moves_prob_vector = moves_vector/sum(moves_vector) ;
moves_prob_vector = cumsum(moves_prob_vector) ; %Tower selector
disp_move_length=13;

%system characteristics
size_of_system=22; %System is considered square shaped
no_of_particles=50;

volume = size_of_system^2;
density = no_of_particles*pi*1*1/volume;


lj_paramA =0.1;
lj_paramB =0.3;

Potential = @ (distance) (lj_paramA/distance^6+lj_paramB/distance^12);

%ProbFunct for GC ensemble
kb=1.38064852e-23;
Prob_Funct = @(DeltaPotentialEnergy,Temperature,Delta_No_particles,No_particles,DNK,mu)...
    ((DNK^Delta_No_particles)*((exp((-DeltaPotentialEnergy+mu*Delta_No_particles)...
    /(kb*Temperature)))*((1/No_particles^Delta_No_particles)+0^(Delta_No_particles>=0)))*No_particles/volume);

temperature = 300;% in kelvin
mu = 0.001;
dnk = .01 ;

kb=1.38064852e-23;%Boltzmann constant

%ProbFunct for Canonical ensemble
%ProbFunct=@(DeltaPotentialEnergy,Temperature)(exp(-DeltaPotentialEnergy/kb*Temperature));
%ProbabilityMove=exp(-DeltaPotentialEnergy/kb*Temperature);

% MicroConfigurations Matrices are defined and initialized

micro_config_matrix = struct('NO_OF_PARTICLES',no_of_particles,'MolecularMatrix',[],'bin_no',[]);

% microstate initialization
micro_config_matrix.MolecularMatrix=[linspace(1,size_of_system-1,no_of_particles);...
    linspace(1,size_of_system-1,no_of_particles)];    

micro_config_matrix.MolecularMatrix=transpose(micro_config_matrix.MolecularMatrix);

sizeMolecularMatrix = size(micro_config_matrix.MolecularMatrix);

% samples are collected
%SampleMemory, will remove it later, only for visualization use in initial stage of the project
sample(attempted_moves) = struct('NO_OF_PARTICLES',[],'sample',zeros(sizeMolecularMatrix),'bin_no',[]) ; 


%Density Measurements and matrix creation

%density mesh defined
density_bin_no = 10 ;
Density=zeros(attempted_moves,density_bin_no,density_bin_no);
density_bin_no = density_bin_no + 1 ;
Density_mesh_x = linspace(0 , size_of_system , density_bin_no) ;
Density_mesh_y = linspace(0 , size_of_system , density_bin_no) ;



for k=1:no_of_particles
    for l=2:density_bin_no
        if micro_config_matrix.MolecularMatrix(k,1) <= Density_mesh_x(l)
            bin_x = l-1;
            break ;
        end
    end
    for m=2:density_bin_no
        if micro_config_matrix.MolecularMatrix(k,2) <= Density_mesh_y(m)
            bin_y = m-1;
            break    ;
        end
    end
    micro_config_matrix.bin_no(k,:) = [bin_x, bin_y ] ;
    Density(1,bin_x ,bin_y) = Density(1,bin_x ,bin_y)+1 ;
end

%here, I naively choose density to visualize convergence
samplesDensity=zeros(attempted_moves,1); 
AverageDensity=zeros(attempted_moves,1);
SumOfDensities=zeros(attempted_moves+1,1);
selected_move_length=0;

%initial density bins defined
micro_config_matrix.NO_OF_PARTICLES = no_of_particles;

% Simulation is started

%transition matrices
transitions_matrix = zeros(200,200);
transition_prob_matrix = zeros(200,200);
bias = 1 ;

% a convergence test
average_energy = zeros(1,attempted_moves);
config_energy = zeros(1,attempted_moves);

for i=1:attempted_moves   
    
    % End Program if the no of particles becomes nil
    if micro_config_matrix.NO_OF_PARTICLES == 0
         input('no more particles')
         return;
    end
    
    % move  selection
     move_selected=zeros(size(moves_prob_vector));
     moves_selector=rand();
     no_move_type = size(moves_prob_vector);
     no_move_type = no_move_type(2);
     
     for k=1:no_move_type
       if moves_selector<=moves_prob_vector(k)
            move_selected(k)=true;
            break;
        end
     end
     
       % Particle Selector
        select=randi(micro_config_matrix.NO_OF_PARTICLES);
        
         % Temporary Move Container is created
    
        temp_config = struct('selectedConf',micro_config_matrix.MolecularMatrix(select,:),...
           'NO_OF_PARTICLES',micro_config_matrix.NO_OF_PARTICLES,'added',false,'removed',false) ;
       
   
    % Now, We make the move
     delta_no_particles = 0;
    
    if (move_selected(1)) %displace Particles
        temp_config.NO_OF_PARTICLES = micro_config_matrix.NO_OF_PARTICLES ;
        
        angle=2*pi*rand(); %move direction select
        selected_move_length=disp_move_length*rand();%move length select
        disp_vector = [selected_move_length*cos(angle),disp_move_length*sin(angle)];
        sizeMolecularMatrix=size(micro_config_matrix);
        
        % temporary variable is updated
        temp_config.selectedConf= temp_config.selectedConf + disp_vector; %movement
        
        % periodic boundary condition
        temp_config.selectedConf=mod(temp_config.selectedConf , [size_of_system,size_of_system]);
        attempted_disp = attempted_disp+1;
    
    elseif (move_selected(2)) %remove Particle
        temp_config.NO_OF_PARTICLES = micro_config_matrix.NO_OF_PARTICLES - 1;
        temp_config.removed = true;
        attempted_remove = attempted_remove+1;
        delta_no_particles = -1;
   
    elseif (move_selected(3)) %add Particle
        
        temp_config.NO_OF_PARTICLES = temp_config.NO_OF_PARTICLES + 1;
        select = temp_config.NO_OF_PARTICLES;
        temp_config.selectedConf=[rand()*size_of_system ,rand()*size_of_system];
        temp_config.added = true;
        attempted_add = attempted_add + 1;
        delta_no_particles = 1;
        
    end
    
    % acceptation criteria
    % Delta Energy Calculation
    delta_potential_energy = 0;
    distance_new = zeros(2);
    distance_old = zeros(2);
   for j = 1 : micro_config_matrix.NO_OF_PARTICLES  % yet unoptimized energy calculation loop
        if j~=select
            if temp_config.added
                distance_new = norm(min(abs(temp_config.selectedConf-micro_config_matrix.MolecularMatrix(j,:)),...
                    abs([size_of_system,size_of_system]-abs(temp_config.selectedConf-micro_config_matrix.MolecularMatrix(j,:)))));
                %new distance of selected particle to j th particle
                distance_old = inf;
                
            elseif temp_config.removed
                distance_new = inf   ;
                distance_old = norm(min(abs(micro_config_matrix.MolecularMatrix(select,:)-micro_config_matrix.MolecularMatrix(j,:))...
                    ,abs([size_of_system,size_of_system]-abs(micro_config_matrix.MolecularMatrix(select,:)-micro_config_matrix.MolecularMatrix(j,:)))));
                %old distance of selected particle to j th particle
                 
                
            else
                distance_new = norm(min(abs(temp_config.selectedConf-micro_config_matrix.MolecularMatrix(j,:)),...
                    abs([size_of_system,size_of_system]-abs(temp_config.selectedConf-micro_config_matrix.MolecularMatrix(j,:)))));
                %new distance of selected particle to j th particle
                distance_old = norm(min(abs(micro_config_matrix.MolecularMatrix(select,:)-micro_config_matrix.MolecularMatrix(j,:))...
                    ,abs([size_of_system,size_of_system]-abs(micro_config_matrix.MolecularMatrix(select,:)-micro_config_matrix.MolecularMatrix(j,:)))));
                %old distance of selected particle to j th particle
            end
            LJpotential = Potential(distance_new) - Potential(distance_old);
            %change in energy by the move for a particle
            delta_potential_energy=delta_potential_energy+LJpotential;
            % net energy difference of the two samples
        end
    end
    
    % (DeltaPotentialEnergy,Temperature,Delta_No_particles,No_particles,DNK,mu)
    probability_move = Prob_Funct(delta_potential_energy,temperature,delta_no_particles,temp_config.NO_OF_PARTICLES,dnk,mu);%for GC Ensemble
    acceptance_prob=min(1,probability_move);
    accept_rand=rand();
   
    %transition matrices updation
     if temp_config.removed  || temp_config.added
            transitions_matrix(micro_config_matrix.NO_OF_PARTICLES,temp_config.NO_OF_PARTICLES) = ...
                transitions_matrix(micro_config_matrix.NO_OF_PARTICLES,temp_config.NO_OF_PARTICLES) + acceptance_prob ;
           
            transitions_matrix(micro_config_matrix.NO_OF_PARTICLES,micro_config_matrix.NO_OF_PARTICLES) = ...
                transitions_matrix(micro_config_matrix.NO_OF_PARTICLES,micro_config_matrix.NO_OF_PARTICLES) + 1 - acceptance_prob ;
     end
    
     if mod(i,20)==0 
        avg_factor = sum(transpose(transitions_matrix));
         for k = 1:200
           transition_prob_matrix(k,:)= transitions_matrix(k,:)/avg_factor(k)^(avg_factor(k)~=0);
%            transition_prob_matrix(isnan(transition_prob_matrix)) = 0;
         end
     end
      bias = transition_prob_matrix(temp_config.NO_OF_PARTICLES,micro_config_matrix.NO_OF_PARTICLES)...transition_prob_matrix(micro_config_matrix.NO_OF_PARTICLES,temp_config.NO_OF_PARTICLES
             /transition_prob_matrix(micro_config_matrix.NO_OF_PARTICLES,temp_config.NO_OF_PARTICLES);
     %bias = 1;
     if isnan(bias)
         bias=1;
     end
      biased_acceptance_prob = bias*acceptance_prob ;
     if i>1
        average_energy(i)= average_energy(i-1);
        config_energy(i)= config_energy(i-1);
     end
     % accepting the move according to Boltzmann statistics
    if accept_rand < biased_acceptance_prob 
        
        config_energy(i) =  config_energy(i) + delta_potential_energy ;
        micro_config_matrix.NO_OF_PARTICLES = temp_config.NO_OF_PARTICLES;
        
        if temp_config.removed
            old_bin = micro_config_matrix.bin_no(select,:) ;
            Density(i,old_bin(1),old_bin(2)) = Density(i,old_bin(1),old_bin(2))- 1 ;
            micro_config_matrix.bin_no(select,:)=[];
            accepted_remove = accepted_remove + 1;
            
        elseif temp_config.added
           micro_config_matrix.MolecularMatrix(select,:) = temp_config.selectedConf;
           for l=2:density_bin_no
                if micro_config_matrix.MolecularMatrix(select,1) <= Density_mesh_x(l)
                    bin_x = l-1;
                    break ;
                end
            end
            for m=2:density_bin_no
                if micro_config_matrix.MolecularMatrix(select,2) <= Density_mesh_y(m)
                    bin_y = m-1;
                    break ;
                end
            end
            micro_config_matrix.bin_no(select,:)= [bin_x, bin_y ] ;
            Density(i,bin_x ,bin_y) = Density(i,bin_x ,bin_y)+1 ;
            accepted_add = accepted_add + 1;
       
         else
            %Bin Change
             old_bin = micro_config_matrix.bin_no(select,:) ;
             Density(i,old_bin(1),old_bin(2)) = Density(i,old_bin(1),old_bin(2))- 1 ;
             micro_config_matrix.MolecularMatrix(select,:) = temp_config.selectedConf;
             for l=2:density_bin_no
                if micro_config_matrix.MolecularMatrix(select,1) <= Density_mesh_x(l)
                    bin_x = l-1;
                    break ;
                end
             end
             for m=2:density_bin_no
                if micro_config_matrix.MolecularMatrix(select,2) <= Density_mesh_y(m)
                    bin_y = m-1;
                    break    ;
                end
             end
            
             micro_config_matrix.bin_no(select,:)= [bin_x, bin_y ] ;
             Density(i,bin_x ,bin_y) = Density(i,bin_x ,bin_y)+1 ;
             accepted_disp = accepted_disp + 1;
           
        end
        accepted_moves=accepted_moves + 1 ;
    end
    average_energy(i) = ((i-1)*average_energy(i) +  config_energy(i))/i;
    
    sample(i).NO_OF_PARTICLES = micro_config_matrix.NO_OF_PARTICLES;% logging, will remove it later
    sample(i).sample = micro_config_matrix.MolecularMatrix;
    sample(i).bin_no = micro_config_matrix.bin_no;
    
    if i<=attempted_moves
        Density(i+1, : , : ) = Density(i, : , : );
    end
    
    % Now, we use density to see the convergence
    no=0; %no of particles in an x X x box
    for k=1:micro_config_matrix.NO_OF_PARTICLES
        if sample(i).sample(k,1)<80 && sample(i).sample(k,2)<80
            no=no+1;
        end
    end
    samplesDensity(i)=no*pi*1*1/(80^2);
    SumOfDensities(i+1)=SumOfDensities(i)+samplesDensity(i);
    AverageDensity(i)=SumOfDensities(i+1)/i;
end

axis equal
scatter(1:attempted_moves,samplesDensity)
plot(1:attempted_moves,AverageDensity)
figure(2)
scatter(sample(100).sample(:,1),sample(100).sample(:,2));
figure(3)
scatter(sample(1000).sample(:,1),sample(1000).sample(:,2));
figure(4)
scatter(sample(5000).sample(:,1),sample(5000).sample(:,2))
figure(5)
scatter(sample(6000).sample(:,1),sample(6000).sample(:,2))
figure(6)
scatter(sample(7000).sample(:,1),sample(7000).sample(:,2))
figure(7)
scatter(sample(8000).sample(:,1),sample(8000).sample(:,2))
figure(8)
scatter(sample(9000).sample(:,1),sample(9000).sample(:,2))
figure(9)
scatter(sample(10000).sample(:,1),sample(10000).sample(:,2))
macro_prob = ones(1,200);
for i=1:199
    macro_prob(i)=  transition_prob_matrix(i,i+1)^(transition_prob_matrix(i,i+1)~=0)/(transition_prob_matrix(i+1,i)^(transition_prob_matrix(i+1,i)~=0));
end
% macro_prob(isnan(macro_prob)) = 1 ;
macro_prob = cumprod(macro_prob);
figure(10)
plot(1:200,macro_prob)
figure(11)
plot(1:attempted_moves,average_energy)
figure(12)
plot(1:attempted_moves,config_energy)


figure(13)
histo_particles=zeros(1,attempted_moves);

for i=1:attempted_moves
    histo_particles(i)=sample(i).NO_OF_PARTICLES;
end
Particles_histogram = histcounts(histo_particles,0:1:200);
plot(1:1:200,Particles_histogram)

    
    
            
                
                
                