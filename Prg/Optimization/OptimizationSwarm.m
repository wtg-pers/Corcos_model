%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Swarm optimization algorythm for MAVLAB
% by N.PICHON 22/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optim_results]=OptimizationSwarm(data,geom,test,obs,method)

e    = 1; %stop conditition

%%%%%Swarm initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n','initialisation...');

CP_C = zeros(test.population,geom.NCP);
CP_T = CP_C;

parfor j1 = 1 : test.population % random initial postion generation
    CP_T(j1,:)      = randi([0 45],1,geom.NCP);
    CP_C(j1,:)      = 0.005+ (0.05-0.005)*rand(1,geom.NCP);
end

parfor j3 = 1 : test.population % random initial speed distribution generation
    V_T(j3,:)       = 20*rand(1,geom.NCP)-10*ones(1,geom.NCP);
    V_C(j3,:)       = 0.02*rand(1,geom.NCP)-0.01*ones(1,geom.NCP);
    
    T_best(j3,:)    = CP_T(j3,:); % record init position as best position know so far to the particle
    C_best(j3,:)    = CP_C(j3,:);
end

parfor j5 = 1 : test.population %calculates aerodynamic/aeroacoustic performance of blade and record it as best so far
    LpBest(j5)      = OptimizationMonoObjective(CP_T(j5,:),CP_C(j5,:),data,geom,test,obs,method);
end

if test.obj(1) == 0 || test.obj(1) == 2 % gets global optimum
    [LpM,Mpos]      = min(LpBest); 
else
    [LpM,Mpos]      = max(LpBest);
end

TM(1:geom.NCP)   = T_best(Mpos,1:geom.NCP); % records absolute best position
CM(1:geom.NCP)   = C_best(Mpos,1:geom.NCP);

parfor j6=1:test.population % calculates particle "distance" in the parameter space to optimum
    dist(j6)        = sqrt(norm(CP_C(j6,:)-CM)^2+10*norm(CP_T(j6,:)-TM)^2);
end
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%work loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LpA       = zeros(1,test.population);
it        = 1;

while max(dist)>e % works while at least one particle is too far from optimum (non-converged)
    fprintf('\n%s\n',['iteration n°',num2str(it)]); %display iteration number
    
    parfor j7 = 1 : test.population
        %changes speed towards optimum position with a random inertia
        r(j7)       = rand;
        V_T(j7,:)   = r(j7)*V_T(j7,:)+(1-r(j7))*(TM-CP_T(j7,:))/(2/0.97725);
        V_C(j7,:)   = r(j7)*V_C(j7,:)+(1-r(j7))*(CM-CP_C(j7,:))/(2/0.97725);
        
        CP_T(j7,:)  = CP_T(j7,:)+V_T(j7,:); %moves particles according to their speed
        CP_C(j7,:)  = CP_C(j7,:)+V_C(j7,:);
    end
    
    for j8 = 1 : test.population
        for i = 1 : geom.NCP %makes sure particles stay within parameter space boundaries
            CP_T(j8,i)  = min(max(0,CP_T(j8,i)),45);
            CP_C(j8,i)  = min(max(0.005,CP_C(j8,i)),0.05);
        end
    end
    
    parfor j9=1:test.population %computes particle performance and record them.        
        LpA(j9)     = OptimizationMonoObjective(CP_T(j9,:),CP_C(j9,:),data,geom,test,obs,method);
        b           = [LpBest(j9),LpA(j9)];
        % Determines best position so far
        if test.obj(1) == 0 || test.obj(1) == 2
            [LpBest(j9),m(j9)]  = min(b);
        else
            [LpBest(j9),m(j9)]  = max(b);
        end
    end
    
    parfor j10=1:test.population %if best position has changed, updates the records      
        if m(j10)==2
            T_best(j10,:) = CP_T(j10,:);
            C_best(j10,:) = CP_C(j10,:);
        end
    end
    
    if test.obj(1) == 0 || test.obj(1) == 2 % gets global optimum
        [LpM,Mpos]  = min(LpBest); 
    else
        [LpM,Mpos]  = max(LpBest);
    end
    
    TM(1:geom.NCP)  = T_best(Mpos,1:geom.NCP); % records best position so far
    CM(1:geom.NCP)  = C_best(Mpos,1:geom.NCP);
    
    parfor j12=1:test.population % calculates particle "distance" in the parameter space to optimum
        dist(j12) = sqrt(norm(CP_C(j12,:)-CM)^2+10*norm(CP_T(j12,:)-TM)^2);
    end
    
    clc;
    it=it+1; %update iteration number
    
end
fprintf(['convergence atteinte en ',num2str(it-1),' itérations']); %display total number of iterations to convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data export%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optim_results.objective_value   = LpM; % records best performance acording to objective
optim_results.geom.T            = TM;  % records best particle geometry
optim_results.geom.C            = CM;
