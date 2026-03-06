%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Multi-objective algorythme based on SPEA 2
%   V2.0 by N.PICHON 22/09/2016
%
%   For more infos see documentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optim_results] = OptimizationSPEA2(data,geom,test,obs,method)
%parpool('AttachedFiles',{'OptimizationMultiObective.m'})
%import OptimizationMultiObjective

%parameters of simulation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N       = test.population;    %size of population (must be multiple of 4)
Na      = N/5;     %size of archive (must be multiple of 4) best practice is N/5
T       = test.generation;     %max number of generations
Ntot    = N+Na ;  %total population

%%init phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%s\n','initialization');

P           = zeros(Ntot,2*geom.NCP);
t           = 1;
K           = floor(sqrt(N+Na));
Pobj        = zeros(Ntot,2);
P_arch      = 0;
Parchold    = 1;
same        = 0;

%random generation of initial population
parfor i1 = 1 : Ntot
    c           = 0.005+0.001 * randi(45,1,geom.NCP);    % chord law
    tw          = randi(45,1,geom.NCP);                  % twist law
    P(i1,:)     = [c tw];                                % initial population in parameter space
end

parfor i10 = 1 : Na
    CP_C         = P(i10,1:geom.NCP);
    CP_T         = P(i10,geom.NCP+1:2*geom.NCP);
    [h,hh]       = OptimizationMultiObjective(CP_T,CP_C,data,geom,test,obs,method);
    Pobj(i10,:)  = [h hh];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%work phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t<=T && same < 3 % works while max number of generations hasn't been reached and pareto front hasn't been stable for 3 generations
    fprintf('\n%s%s\n','iteration n°',num2str(t));
%%%%%projection of population into the objective space
    fprintf('%s\n','projection into objective space');
    parfor i2 = Na+1 : Ntot
        CP_C         = P(i2,1:geom.NCP);
        CP_T         = P(i2,geom.NCP+1:2*geom.NCP);
        %addAttachedFiles((Ntot-Na+1), 'OptimizationMultiObective.m');
        %addAttachedFiles(myPool, 'myFunction1.m');
        %addAttachedFiles(poolobj,{'myFun1.m','myFun2.m'})
        [a100,b100]  = OptimizationMultiObective(CP_T,CP_C,data,geom,test,obs,method);
        Pobj(i2,:)   = [a100,b100];
    end
%%%%%domination check
    fprintf('%s\n','domination check');
    domination      = zeros (Ntot,Ntot);
    is_dominated    = zeros(1,Ntot);
    %domination relations matrix. domination(i,j) = 1 means i dominates j, 0 means non-dominates, -1 is dominated by
    for i3 = 1 : Ntot-1
        for j = i3+1 : Ntot
            if Pobj(i3,1) <= Pobj(j,1) && Pobj(i3,2) <= Pobj(j,2)
                if Pobj(i3,1) < Pobj(j,1)
                    domination(i3,j)    = 1;
                    is_dominated(j)     = 1;
                elseif Pobj(i3,2) < Pobj(j,2)
                    domination(i3,j)    = 1;
                    is_dominated(j)     = 1;
                else
                    domination(i3,j)    = 0;
                end
            end
            if Pobj(i3,1) >= Pobj(j,1) && Pobj(i3,2) >= Pobj(j,2)
                if Pobj(i3,1) > Pobj(j,1)
                    domination(i3,j)    = -1;
                    is_dominated(i3)    = 1;
                elseif Pobj(i3,2) > Pobj(j,2)
                    domination(i3,j)    = -1;
                    is_dominated(i3)    = 1;
                else
                    domination(i3,j)    = 0;
                end
            end
        end
    end
    domination = domination + domination';
    
%%%%%strength calculation
    fprintf('%s\n','strength and fitness calculation');
    D       = zeros(Ntot,Ntot);
    sigma   = zeros(Ntot);
    
    parfor i4 = 1 : Ntot % Strength calculation
        S(i4)   = 0;
        for j2 = 1 : Ntot
            if domination(i4,j2) == 1
                S(i4)   = S(i4) +1/(Ntot+1);
            end
        end
    end
    for i8 = 1 : Ntot-1 % distance calculation in the objective space
        for j4 = i8 : Ntot
            D(i8,j4)    = sqrt((Pobj(i8,1)-Pobj(j4,1))^2+(Pobj(i8,2)-Pobj(j4,2))^2);
        end
    end
    D   = D+D'; %symetrization of D
    D   = sort(D); %sorting Distance
    parfor i9 = 1 : Ntot
        sigma(i9)       = D(i9,K); %calculation of sigma, the k-th distance vector
    end
    
    
    
    parfor i5 = 1 : Ntot %fitness calculation
        F(i5)   = 0;
        for j3 = 1 : Ntot
            if domination(i5,j3) == 1
                F(i5)   = F(i5) +S(j3)+1/(sigma(j3)+2);
            end
        end
    end
    
    clear Parchold;
    Parchold    = P_arch;
    clear P_arch;
    clear P_arch_obj;
    P_arch      = [];
    P_obj_arch  = [];
%%%%%archive management
    fprintf('%s\n','archive management');
    parfor i6 = 1 : Ntot  %take all non-dominated to archive
        if is_dominated(i6) == 0
            P_arch      = [P_arch ; P(i6,:)];
            P_obj_arch  = [P_obj_arch ; Pobj(i6,:)];
        end
    end
    clear A; clear A_obj;
    A       = P_arch;
    A_obj   = P_obj_arch;
    
    if size(P_arch,1) > Na % if archive too big, truncate
        [P_arch,P_obj_arch] = OptimizationSPEA2Truncate(P_arch,P_obj_arch,D,Na,Ntot);
    end
    
    
    
    if size(P_arch,1) < Na
        k = size(P_arch,1);
        for i7 = 1 : Na-k
            P_arch      = [P_arch ; P(k+i7,:)];
            P_obj_arch  = [P_obj_arch ; Pobj(k+i7,:)];
        end
    end
    
%%%%%mating
    fprintf('%s\n','mating');
    Pmat    = zeros(N/2,2*geom.NCP);
    
    for i=1:N/2 %binary tournament to chose mating pop
        a=randi(length(P));
        b=a;
        while b==a
            b=randi(length(P));
        end
        
        if F(a)>F(b) %selection of the fitest
            Pmat(i,:)=P(a,:);
            P(a,:)=[];
            if b > a
                P(b-1,:) = [];
            else
                P(b,:) = [];
            end
        else
            Pmat(i,:)=P(b,:);
            P(a,:)=[];
            if b > a
                P(b-1,:) = [];
            else
                P(b,:) = [];
            end
        end
    end
    
    Pnew           = zeros(Ntot,2*geom.NCP);
    Pobjnew        = zeros(Ntot,2);
    Pnew(1:Na,:)   = P_arch; %save archive as part of the new population
    Pobjnew(1:Na,:)   = P_obj_arch;
    
    for i100=1:N/4 %mating of random couples taken into the mating pool
        a=randi(size(Pmat,1));
        b=a;
        while b==a
            b=randi(size(Pmat,1));
        end
        %4 children are generated with one point crossover
        c = randi([2 2*geom.NCP-1],1,4);
        for j100=1:4
            Pnew(Na+4*(i100-1)+j100,1:c(j100)) = Pmat(a,1:c(j100));
            Pnew(Na+4*(i100-1)+j100,c(j100)+1:end) = Pmat(b,c(j100)+1:end);
            for k100 = 1 : geom.NCP % include a 5% mutation rate
                r2 = rand;
                if r2>0.95
                    Pnew(Na+4*(i100-1)+j100,k100)= 0.005 + 0.001*randi(45);
                end
            end
            for k100 = geom.NCP+1 : 2*geom.NCP
                r3 = rand;
                if r3>0.95
                    Pnew(Na+4*(i100-1)+j100,k100)= randi(45);
                end
            end
        end
    end
    
    clear P;
    clear Pobj;
    P        = Pnew;
    Pobj     = Pobjnew;
    clear Pnew;
    clear Pobjnew;
%%%%%%%%%%
    t=t+1;
    if isequal(Parchold,P_arch)==1
        same = same +1;
    else
        same = 0;
    end
    clc;
end
%%%%%data export%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optim_results.objective_value   = A_obj; % records best performance acording to objectives
optim_results.geom              = A;     % records best particle geometry

