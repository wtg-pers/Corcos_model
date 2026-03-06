%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User-selected Combinatory optimization
%
% V1.0 by R.SERRE and N.PICHON 22/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optim_results] = OptimizationCombination(data,geom,test,obs,method)

%initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial range for twist and chord
T_max       = 30;
T_min       = 5;

T_range     = T_min : (T_max-T_min)/(test.grid-1) : T_max;

T_set       = {};
for i1 = 1 : geom.NCP-1
    T_set{1,i1}    =  T_range;
end

%initial blade table that is made of all possible combinations
switch geom.NCP
    case 3
        [b1, b2]        = ndgrid(T_set{:});
        CP_T = [b1(:) b2(:)];
    case 4
        [b1, b2, b3]    = ndgrid(T_set{:});
        CP_T = [b1(:) b2(:) b3(:)];
    case 5
        [b1, b2, b3, b4]    = ndgrid(T_set{:});
        CP_T = [b1(:) b2(:) b3(:) b4(:)];
    case 6
        [b1, b2, b3, b4, b5]    = ndgrid(T_set{:});
        CP_T = [b1(:) b2(:) b3(:) b4(:) b5(:)];
    case 7
        [b1, b2, b3, b4, b5, b6]    = ndgrid(T_set{:});
        CP_T = [b1(:) b2(:) b3(:) b4(:) b5(:) b6(:)];
    case 8
        [b1, b2, b3, b4, b5, b6, b7]    = ndgrid(T_set{:});
        CP_T = [b1(:) b2(:) b3(:) b4(:) b5(:) b6(:) b7(:)];
end

CP_T        = [CP_T zeros(size(CP_T,1),1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OBJ = zeros(2,size(CP_T,1));
for j = 1 : size(CP_T,1) %compute efficiency relative to objectives
    fprintf('%s\n',['remaining blades to test : ',num2str(size(CP_T,1)-(j-1)),...
        '/',num2str(size(CP_T,1))]);
    [a,b] = OptimizationMultiObjective(CP_T(j,:),geom.CP_C(:)',data,geom,test,obs,method);
    OBJ(:,j) = [a b];
end

%find pareto front%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n%s\n','computing Pareto Front...');
is_dominated = zeros(size(CP_T,1));
    for j1 = 1 : size(CP_T,1)
            for k2 = 1 : size(CP_T,1)
                if (j1~=k2)
                if OBJ(1,j1)<OBJ(1,k2) && OBJ(2,j1)<=OBJ(2,k2)
                    is_dominated(k2) = 1;
                elseif OBJ(1,j1)<=OBJ(1,k2) && OBJ(2,j1)<OBJ(2,k2)
                    is_dominated(k2) = 1;
                elseif OBJ(1,j1)>=OBJ(1,k2) && OBJ(2,j1)>OBJ(2,k2)
                    is_dominated(j1) = 1;
                elseif OBJ(1,j1)>OBJ(1,k2) && OBJ(2,j1)>=OBJ(2,k2)
                    is_dominated(j1) = 1;
                end
                end
            end
    end
    %for j1 = 1 : size(CP_T,1)-1
    %        for k2 = j1+1 : size(CP_T,1)
    %            if OBJ(1,j1)>OBJ(1,k2) && OBJ(2,j1)>=OBJ(2,k2)
    %                is_dominated(k2) = 1;
    %            elseif OBJ(1,j1)>=OBJ(1,k2) && OBJ(2,j1)>OBJ(2,k2)
    %                is_dominated(k2) = 1;
    %            elseif OBJ(1,j1)<=OBJ(1,k2) && OBJ(2,j1)<OBJ(2,k2)
    %                is_dominated(j1) = 1;
    %            elseif OBJ(1,j1)<OBJ(1,k2) && OBJ(2,j1)<=OBJ(2,k2)
    %                is_dominated(j1) = 1;
    %            end
    %        end
    %end

pareto = 1-is_dominated;

pareto_obj = []; % register performance and position in chord/twist table of pareto front elements
    for j1 = 1 : size(CP_T,1)
        if pareto(j1) == 1
            a = [OBJ(1,j1) OBJ(2,j1) j1];
            pareto_obj = [pareto_obj; a];
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
    optim_results.geom.T(i1,1:geom.NCP)  = CP_T(pareto_obj(i1,3),:);
    %optim_results.geom.T(i1,1:geom.NCP)  = CP_T(pareto_obj(i1,4),:);
end
