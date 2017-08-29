%made by Ishan Arora under Dr.Rane, IITGN, Gandhinagar
%We use 6-12 potentials, 
%optimization tactics not employed yet
%also need to improve readability of the code
%replace by suitable(reduced)units later.
moveLength=5;%maximum move length
HardDiskDia=2;%Dont need it yet
NO_OF_ITERATIONS=100000;
SIZE_OF_SYSTEM=220;%length of square 2D box
NO_OF_PARTICLE=90;
lj_paramA=.1;%lj parameters
lj_paramB=.3;
Temperature=300;% in kelvin
kb=1.38064852e-23;%Boltzmann constant
MolecularMatrix=[linspace(1,SIZE_OF_SYSTEM-1,NO_OF_PARTICLE)...
    ;linspace(1,SIZE_OF_SYSTEM-1,NO_OF_PARTICLE)];%initial Locations of n particles arranged diagonally across the square box
MolecularMatrix=transpose(MolecularMatrix);%transposed to allow n position vectors of particles to be in columns as expected for efficiency 
sizeMolecularMatrix=size(MolecularMatrix);
sample=ones(NO_OF_ITERATIONS,sizeMolecularMatrix(1),sizeMolecularMatrix(2)) ;%SampleMemory, will remove it later, only for visualization use in initial stage of the project
%density average
VOLUME=SIZE_OF_SYSTEM^2;
density=NO_OF_PARTICLE*pi*1*1/VOLUME;
samplesDensity=zeros(NO_OF_ITERATIONS,1); %here, I naively choose density to visualize convergence
AverageDensity=zeros(NO_OF_ITERATIONS,1);
SumOfDensities=zeros(NO_OF_ITERATIONS+1,1);
selectedMoveLength=0;
acceptedMoves=0;
for i=1:NO_OF_ITERATIONS   %sampleGenerator i->no of Samples
    %I need to implement better method for random sampling
     angle=2*pi*rand(); %move direction select
     selectedMoveLength=moveLength*rand();%move length select
     d=[selectedMoveLength*cos(angle),moveLength*sin(angle)];
     sizeMolecularMatrix=size(MolecularMatrix);
     select=randi(NO_OF_PARTICLE,1);           %selection
     BackUpValue=MolecularMatrix(select,:);
     BackUpValue= BackUpValue +d;          %movement
   %....................................PERIODIC BOUNDARY
     BackUpValue=mod(BackUpValue,[SIZE_OF_SYSTEM,SIZE_OF_SYSTEM]);
     DeltaPotentialEnergy=0;
    delta_new=zeros(2);
    delta_old=zeros(2);
    for j=1:NO_OF_PARTICLE%yet unoptimized energy calculation loop
        if j~=select
           delta_new=min(abs(BackUpValue-MolecularMatrix(j,:)),...
                abs([SIZE_OF_SYSTEM,SIZE_OF_SYSTEM]-abs(BackUpValue-MolecularMatrix(j,:))));%new distance of selected particle to j th particle
            delta_old=min(abs(MolecularMatrix(select,:)-MolecularMatrix(j,:))...
                ,abs([SIZE_OF_SYSTEM,SIZE_OF_SYSTEM]-abs(MolecularMatrix(select,:)-MolecularMatrix(j,:))));%old distance of selected particle to j th particle
            LJpotential=(lj_paramA/(norm(delta_new))^6+lj_paramB/(norm(delta_new))^12)...
                -(lj_paramA/(norm(delta_old))^6+lj_paramB/(norm(delta_old))^12);%change in energy by the move for a particle
           DeltaPotentialEnergy=DeltaPotentialEnergy+LJpotential;% net energy difference of the two sampless
        end
    end
    BoltzmannFactor=exp(-DeltaPotentialEnergy/kb*Temperature);
    AcceptanceProb=min(1,BoltzmannFactor);
    selectRand=rand();
    if selectRand<AcceptanceProb %accepting the move according to Boltzmann statistics
        MolecularMatrix(select,:,:)=BackUpValue;   
        acceptedMoves=acceptedMoves+1;
    end
    sample(i,:,:)=MolecularMatrix;% will remove it later
    %.........................................................Now, we use density to see the convergence 
    no=0; %no of particles in an x X x box
    for k=1:NO_OF_PARTICLE
        if sample(i,k,1)<80 && sample(i,k,2)<80
            no=no+1;
        end
    end
    samplesDensity(i)=no*pi*1*1/(80^2);
    SumOfDensities(i+1)=SumOfDensities(i)+samplesDensity(i);
    AverageDensity(i)=SumOfDensities(i+1)/i;
end
subplot(2,1,1)
scatter([1:NO_OF_ITERATIONS],samplesDensity)
subplot(2,1,2)
plot([1:NO_OF_ITERATIONS],AverageDensity)
figure(2)
scatter(sample(1,:,1),sample(1,:,2))
figure(3)
scatter(sample(100,:,1),sample(100,:,2))
figure(4)
scatter(sample(1000,:,1),sample(1000,:,2))
figure(5)
scatter(sample(10000,:,1),sample(10000,:,2))



    
            
                
                
                
