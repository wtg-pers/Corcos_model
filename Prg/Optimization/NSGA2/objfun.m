function [y, cons] = objfun( pop, geom,test,obs,method )

geom.CP_C(1:geom.NCP) = pop(1:geom.NCP);
geom.CP_T(1:geom.NCP) = pop(geom.NCP+1:2*geom.NCP);


[~,~,flow,method,~,~,test,acoustic] = MAVLAB (geom,test,obs,method);

if test.eliminate == 0
    if test.obj(1) == 0
        if method.dBA == 1
            OBJ_value = acoustic.OASPLA_tot;
        else
            OBJ_value = acoustic.OASPL_tot;
        end
    elseif test.obj(1) == 1
        OBJ_value = -flow.Thrust;
    else
        OBJ_value = flow.P;
    end
    
    
    if test.obj(2) == 0
        if method.dBA == 1
            OBJ_value_2 = acoustic.OASPLA_tot;
        else
            OBJ_value_2 = acoustic.OASPL_tot;
        end
    elseif test.obj(2) == 1
        OBJ_value_2 = -flow.Thrust;
    else
        OBJ_value_2 = flow.P;
    end
else
    if test.obj(1) == 0
        if method.dBA == 1
            OBJ_value = 1000000000000000;
        else
            OBJ_value = 1000000000000000;
        end
    elseif test.obj(1) == 1
        OBJ_value = -1000000000000000;
    else
        OBJ_value = 1000000000000000;
    end
    
    
    if test.obj(2) == 0
        if method.dBA == 1
            OBJ_value_2 = 1000000000000000;
        else
            OBJ_value_2 = 1000000000000000;
        end
    elseif test.obj(2) == 1
        OBJ_value_2 = -1000000000000000;
    else
        OBJ_value_2 = 1000000000000000;
    end
end


y(1)  = OBJ_value;
y(2)  = OBJ_value_2;

cons = [];