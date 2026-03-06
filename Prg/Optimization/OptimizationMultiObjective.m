%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiobj_adapter
% V1.0 by N.PICHON 22/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OBJ_value_1 OBJ_value_2] = OptimizationMultiObjective(CP_T,CP_C,data,geom,test,obs,method)

geom.CP_C = CP_C;
geom.CP_T = CP_T;

[~,~,flow,method,~,~,test,acoustic] = RotorEvaluation(data,geom,test,obs,method);

if test.eliminate == 0
    if test.obj(1) == 0
        if method.dBA == 1
            OBJ_value_1 = acoustic.OASPLA_Total;
        else
            OBJ_value_1 = acoustic.OASPL_Total;
        end
    elseif test.obj(1) == 1
        OBJ_value_1 = flow.Thrust;
    else
        OBJ_value_1 = flow.P;
    end
    
    if test.obj(2) == 0
        if method.dBA == 1
            OBJ_value_2 = acoustic.OASPLA_Total;
        else
            OBJ_value_2 = acoustic.OASPL_Total;
        end
    elseif test.obj(2) == 1
        OBJ_value_2 = flow.Thrust;
    else
        OBJ_value_2 = flow.P;
    end
else
    if test.obj(1) == 0
        if method.dBA == 1
            OBJ_value_1 = 100000000000000000;
        else
            OBJ_value_1 = 100000000000000000;
        end
    elseif test.obj(1) == 1
        OBJ_value_1 = 0;
    else
        OBJ_value_1 = 100000000000000000;
    end
    
    if test.obj(2) == 0
        if method.dBA == 1
            OBJ_value_2 = 100000000000000000;
        else
            OBJ_value_2 = 100000000000000000;
        end
    elseif test.obj(2) == 1
        OBJ_value_2 = 0;
    else
        OBJ_value_2 = 100000000000000000;
    end
end
