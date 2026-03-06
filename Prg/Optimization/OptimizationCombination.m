%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-selected Combinatory optimization
%
% V1.0 by R.SERRE and N.PICHON 22/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optim_results] = OptimizationCombination(data,geom,test,obs,method)

%initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial range for twist and chord
T_max       = test.T_max;
T_min       = test.T_min;
C_max       = test.C_max;
C_min       = test.C_min;

C_range     = C_min : (C_max-C_min)/(test.grid-1) : C_max;
T_range     = T_min : (T_max-T_min)/(test.grid-1) : T_max;

C_set       = {};
T_set       = {};
for i1 = 1 : geom.NCP
    C_set{1,i1}    =  C_range;
end
for i1 = 1 : geom.NCP-1
    T_set{1,i1}    =  T_range;
end

%initial blade table that is made of all possible combinations
switch geom.NCP
    case 3
        [a1, a2, a3]        = ndgrid(C_set{:});
        [b1, b2]        = ndgrid(T_set{:});
        CP_C = [a1(:) a2(:) a3(:)];
        CP_T = [b1(:) b2(:)];
    case 4
        [a1, a2, a3, a4]    = ndgrid(C_set{:});
        [b1, b2, b3]    = ndgrid(T_set{:});
        CP_C = [a1(:) a2(:) a3(:) a4(:)];
        CP_T = [b1(:) b2(:) b3(:)];
    case 5
        [a1, a2, a3, a4, a5]    = ndgrid(C_set{:});
        [b1, b2, b3, b4]    = ndgrid(T_set{:});
        CP_C = [a1(:) a2(:) a3(:) a4(:) a5(:)];
        CP_T = [b1(:) b2(:) b3(:) b4(:)];
    case 6
        [a1, a2, a3, a4, a5, a6]    = ndgrid(C_set{:});
        [b1, b2, b3, b4, b5]    = ndgrid(T_set{:});
        CP_C = [a1(:) a2(:) a3(:) a4(:) a5(:) a6(:)];
        CP_T = [b1(:) b2(:) b3(:) b4(:) b5(:)];
    case 7
        [a1, a2, a3, a4, a5, a6, a7]    = ndgrid(C_set{:});
        [b1, b2, b3, b4, b5, b6]    = ndgrid(T_set{:});
        CP_C = [a1(:) a2(:) a3(:) a4(:) a5(:) a6(:) a7(:)];
        CP_T = [b1(:) b2(:) b3(:) b4(:) b5(:) b6(:)];
    case 8
        [a1, a2, a3, a4, a5, a6, a7, a8]    = ndgrid(C_set{:});
        [b1, b2, b3, b4, b5, b6, b7]    = ndgrid(T_set{:});
        CP_C = [a1(:) a2(:) a3(:) a4(:) a5(:) a6(:) a7(:) a8(:)];
        CP_T = [b1(:) b2(:) b3(:) b4(:) b5(:) b6(:) b7(:)];
end

CP_T        = [CP_T zeros(size(CP_T,1),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
p = parpool(test.CoreNumber);
%
OBJ = zeros(2,size(CP_C,1),size(CP_T,1));
%
ilocal = 1 ; 
if(exist('./BackUpFolder/BackUp.mat','file')>0)
load './BackUpFolder/BackUp.mat';
end
% the variable ilocal was saved and reset to continue the calculations from where it was left off
for i = ilocal : size(CP_C,1)
    fprintf('%s\n',['remaining blades to test : ',num2str(size(CP_T,1)*size(CP_C,1)-(i-1)*size(CP_T,1)),...
        '/',num2str(size(CP_T,1)*size(CP_C,1))]);
    parfor j = 1 : size(CP_T,1) %compute efficiency relative to objectives
        [a,b] = OptimizationMultiObjective(CP_T(j,:),CP_C(i,:),data,geom,test,obs,method);
        OBJ(:,i,j) = [a b];
    end
    ilocal=i;
    save('./BackUpFolder/BackUp.mat','ilocal','OBJ');
end
%p.NumWorkers
%
%MAX=1000000;
%i=1;
%if(exist('result.mat','file')>0)
%load 'result.mat';
%end
%% the variable i was saved and reset to continue the calculations from
%where it was left off
%for n=i:MAX
%docalculations;
%save('result.mat','i','variable1','variable2',..);
%end
%
%find pareto front%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n%s\n','computing Pareto Front...');
is_dominated = zeros(size(CP_C,1),size(CP_T,1));
for i1 = 1 : size(CP_C,1)-1
    for j1 = 1 : size(CP_T,1)-1
        for k1 = i1+1 : size(CP_C,1)
            for k2 = j1+1 : size(CP_T,1)
                if OBJ(1,i1,j1)<OBJ(1,k1,k2) && OBJ(2,i1,j1)<=OBJ(2,k1,k2)
                    is_dominated(k1,k2) = 1;
                elseif OBJ(1,i1,j1)<=OBJ(1,k1,k2) && OBJ(2,i1,j1)<OBJ(2,k1,k2)
                    is_dominated(k1,k2) = 1;
                elseif OBJ(1,i1,j1)>=OBJ(1,k1,k2) && OBJ(2,i1,j1)>OBJ(2,k1,k2)
                    is_dominated(i1,j1) = 1;
                elseif OBJ(1,i1,j1)>OBJ(1,k1,k2) && OBJ(2,i1,j1)>=OBJ(2,k1,k2)
                    is_dominated(i1,j1) = 1;
                end
            end
        end
    end
end

pareto = 1-is_dominated;

pareto_obj = []; % register performance and position in chord/twist table of pareto front elements
for i1 = 1 : size(CP_C,1)
    for j1 = 1 : size(CP_T,1)
        if pareto(i1,j1) == 1
            a = [OBJ(1,i1,j1) OBJ(2,i1,j1) i1 j1];
            pareto_obj = [pareto_obj; a];
        end
    end
end


%%plot pareto_front%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure; % plot pareto front efficiency
%plot(pareto_obj(:,2),pareto_obj(:,1),'.','linewidth',2);
%hold on;
%
%switch test.obj(1); % label axes in accordance with objectives
%    case 0
%        ylabel('Overall Noise');
%    case 1
%        ylabel('Thrust');
%    case 2
%        ylabel('Power Consumption');
%end
%switch test.obj(2);
%    case 0
%        xlabel('Overall Noise');
%    case 1
%        xlabel('Thrust');
%    case 2
%        xlabel('Power Consumption');
%end

%export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optim_results.objective_value        = pareto_obj(:,1:2);
for i1 = 1 : size(pareto_obj,1)
    optim_results.geom.C(i1,1:geom.NCP)  = CP_C(pareto_obj(i1,3),:);
    optim_results.geom.T(i1,1:geom.NCP)  = CP_T(pareto_obj(i1,4),:);
end
