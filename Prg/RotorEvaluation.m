%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAVLAB                                                                  %
%                                                                         %
% Framework for aeroacoustic optimization tool                            %
%                                                                         %
%  Main contributors ========                                             %
%  C. NANA          (2015)                                                %
%  A. LEMAISTRE     (2015)                                                %
%  T. JARDIN        (2015)                                                %
%  L. CAMPO-CARMONA (2016)                                                %
%  N. PICHON        (2016)                                                %
%  V. CHAPIN        (2016)                                                %
%  R. SERRE         (2017)                                                %
%  Latest contributor =======                                             %
%  R. SERRE (21/03/2017)                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [geom,data,flow,method,obs,polar,test,acoustic] = RotorEvaluation(data,geom,test,obs,method)
%
%addpath(genpath('./Prg/'));
%
if method.optim == -1
    fprintf('\n%s\n','Initialization');
end
%
% GEOMETRY ==================================================================================
if (method.custom==1) %custom geometry definition
    %
    dCPR      = (geom.Rtip-geom.Rhub)/(geom.NCP-1); % radial distance between two CP
    geom.CP_R = zeros(1,geom.NCP);
    %
    if (method.airfopt==0)
    geom.CP_R   = geom.Rhub:dCPR:geom.Rtip;         %radius of the NCP control points, equally distributed on the blade
    elseif (method.airfopt==1)
    geom.CP_R(1)   = geom.Rhub;
    geom.CP_R(end) = geom.Rtip;
    dCP_airfoil = (geom.CP_R(end)-geom.CP_R(1))/(method.NFoilOpt-1);CP_airfoil=zeros(method.NFoilOpt,1);
    CP_airfoil  = [geom.CP_R(1) : dCP_airfoil : geom.CP_R(end)];
    % ensure that the radial control points match airfoil section discretisation
    for CP=2:geom.NCP-1
       geom.CP_R(CP) = CP_airfoil( find( abs(CP_airfoil(:)-(CP_airfoil(1)+(CP-1)*dCPR))==min(abs(CP_airfoil(:)-(CP_airfoil(1)+(CP-1)*dCPR))) ) );
    end
    end 
    dRq         = (geom.Rtip-geom.Rhub)/(geom.Nq-1);  % radial distance between two elements centers
    geom.Rq     = geom.Rhub:dRq:geom.Rtip;          %radius of blade elements centers
    geom.dR     = dRq;
    if (method.sweep==0)
        geom.CP_Fx = zeros(1,geom.NCP); %Sweep X of CP
    end
    switch method.interpolation %interpolation of blade elements centers chord, twist and sweep X
        case 0
            geom.C  = interp1(geom.CP_R,geom.CP_C,geom.Rq,'spline');
            geom.T  = interp1(geom.CP_R,geom.CP_T,geom.Rq,'spline');
            geom.Fx = interp1(geom.CP_R,geom.CP_Fx,geom.Rq,'spline');
        case 1
            geom.C  = Geom_bezier(geom.CP_R,geom.CP_C,geom.Rq);
            geom.T  = Geom_bezier(geom.CP_R,geom.CP_T,geom.Rq);
            geom.Fx = Geom_bezier(geom.CP_R,geom.CP_Fx,geom.Rq);
    end
elseif (method.custom==0) %reference geometry loading
    fid = fopen(['./geometry_files/',geom.blade_ref,'.geom'],'r'); % loading of reference geometry file
    M = textscan(fid,'%f %f %f %f %f','Delimiter','\t','Headerlines',4); % file scan
    fclose(fid);
    geom.Rq = M{1}; %radius of blade elements centers
    geom.C  = M{2}; %chord of blade elements centers
    geom.T  = M{3}; %twist of blade elements centers
    
    if (method.sweep==1)
        geom.Fx = M{4}; %sweep X of blade elements centers
    else
        geom.Fx = zeros(1,length(geom.Rq)); %sweep X of blade elements centers
    end
    geom.Nq     = length(geom.Rq); %number of blade elements
    geom.Rhub   = geom.Rq(1); %radius of blade hub
    geom.Rtip   = geom.Rq(end); %radius of blade hub
    geom.dR = geom.Rq(2)-geom.Rq(1);
elseif (method.custom==2)
    FilterIndex = 0;
    while FilterIndex == 0
    [FileName,PathName,FilterIndex] = uigetfile('*.geom','select geometry file');
    end
    fid = fopen([PathName,FileName],'r'); % loading of reference geometry file
    M = textscan(fid,'%f %f %f %f %f','Delimiter','\t','Headerlines',4); % file scan
    fclose(fid);
    geom.Rq = M{1}; %radius of blade elements centers
    geom.C  = M{2}; %chord of blade elements centers
    geom.T  = M{3}; %twist of blade elements centers
    
    if (method.sweep==1)
        geom.Fx = M{4}; %sweep X of blade elements centers
    else
        geom.Fx = zeros(1,length(geom.Rq)); %sweep X of blade elements centers
    end
    geom.Nq     = length(geom.Rq); %number of blade elements
    geom.Rhub   = geom.Rq(1); %radius of blade hub
    geom.Rtip   = geom.Rq(end); %radius of blade hub
    geom.dR = geom.Rq(2)-geom.Rq(1);
end
%
if (method.weight==1) % Add the weight of the blade to the total mass to lift
    if (method.airfopt==0)
    geom.Weight = Aerodynamic_MonoAirfoil_WeightWatcher(geom.airfoil,geom.Rq,geom.C,geom.NBlades);
    elseif (method.airfopt==1)
    geom.Weight = Aerodynamic_MultiAirfoil_WeightWatcher(geom.Rq,geom.C,geom.NBlades,method);
    end 
    test.thrust  = test.thrust + geom.Weight* 9.81;
end
%
% AERODYNAMIC ===============================================================================
% polar
if (method.polar==1)
    if (method.airfopt==0)
    polar = Aerodynamic_MonoAirfoil_PolarSurface(geom,method); %polar surface for future interpolation creation
    elseif (method.airfopt==1)
    polar = Aerodynamic_MultiAirfoil_PolarSurface(geom,method); %polar surface for future interpolation creation
    end 
elseif (method.polar==0)
    polar_file = ['./PolarFunction/',geom.airfoil,'.fun'];
    polar = Aerodynamic_polarfunction(polar_file);
end
if (method.optim==-1)
    fprintf('%s\n','Aerodynamic computation');
end
% BEMT
if (method.isothrust==1);
    if (method.airfopt==0)
    [flow,test,geom] = Aerodynamic_MonoAirfoil_IsoThrustV1(test,geom,polar,data,method); %iso-thrust test at specified thrust goal
    elseif (method.airfopt==1)
    [flow,test,geom] = Aerodynamic_MultiAirfoil_IsoThrustV1(test,geom,polar,data,method); %iso-thrust test at specified thrust goal
    end 
elseif (method.isothrust==2);
    if (method.airfopt==0)
    [flow,test,geom] = Aerodynamic_MonoAirfoil_IsoThrustV2(test,geom,polar,data,method); %iso-thrust test at specified thrust goal
    elseif (method.airfopt==1)
    [flow,test,geom] = Aerodynamic_MultiAirfoil_IsoThrustV2(test,geom,polar,data,method); %iso-thrust test at specified thrust goal
    end 
else (method.isothrust==0);
%fprintf('Hello\n')
    if (method.airfopt==0)
    [flow,geom]      = Aerodynamic_MonoAirfoil_BEMT(test,geom,polar,data,method); %BEMT test at specified RPM
    elseif (method.airfopt==1)
    [flow,geom]      = Aerodynamic_MultiAirfoil_BEMT(test,geom,polar,data,method); %BEMT test at specified RPM
    end 
end
%
geom.D   = 2*max(geom.Rb);
n        = test.RPM/60;
flow.P   = 2*pi*n*flow.Torque;

flow.CT   = flow.Thrust/(data.rho* (pi*(geom.D/2)^2) * ((n*2*pi) * geom.D/2)^2)
flow.CP   = flow.P/(data.rho* (pi*(geom.D/2)^2) * ((n*2*pi) * geom.D/2)^3)
flow.CQ   = flow.CP/(n*2*pi)
%flow.CT  = flow.Thrust/(data.rho*n^2*geom.D^4);
%flow.CP  = flow.P /(data.rho*n^3*geom.D^5);
%flow.CQ  = flow.CP/(2*pi);
%test blade validity
if (method.blade_validity==1)
    [test.eliminate]     = Aerodynamic_DetachmentDetector(geom,flow,test);
else
    test.eliminate       = 0;
end
%
% ACOUSTIC ==================================================================================
% OPTIMIZATION ==============================================================================
acoustic = [];
if (method.acoustic==1)
    if (test.eliminate==0)
        if (test.error==0)
            if (method.optim==-1)
                fprintf('%s\n','Tonal noise computation');
            end
            [obs,acoustic]      = Acoustic_tonal(test,geom,flow,obs,data,method); %tonal noise through Casalino's resolution of Ffowcs-Williams and Hawking Equation
            acoustic.OASPL_Total  = acoustic.OASPL_Tonal;
            acoustic.Lp_Total     = acoustic.Lp_Tonal';
            acoustic.OASPLA_Total = acoustic.OASPLA_Tonal;
            acoustic.LpA_Total    = acoustic.LpA_Tonal';
            if (method.freqband==0)
               acoustic.fbb = acoustic.f;
            elseif (method.freqband==1)
               acoustic.fbb = acoustic.foct;
            end
            if (method.broadband==1)
                if (method.optim==-1)
                    fprintf('%s\n','Broadband noise computation (trailing-edge)');
                end
                acoustic            = Acoustic_broadband_TE(test,geom,flow,acoustic,obs,data,method); %broadband noise through Rozenberg et al's method
                %
                if (method.freqband==0)
                   acoustic.OASPL_tot    = 10*log10(10^(0.1*acoustic.OASPL_Total)+10^(0.1*acoustic.OASPL_BB_TE));
                   acoustic.Lp_Total     = 10*log10(10.^(0.1*acoustic.Spp_BB_TE') + 10.^(0.1*acoustic.Lp_Total));
                elseif (method.freqband==1)
                   acoustic.OASPLA_Total = 10*log10(10^(0.1*acoustic.OASPLA_Total)+10^(0.1*acoustic.OASPLA_BB_TE));
                   acoustic.LpA_Total    = 10*log10(10.^(0.1*acoustic.SppA_BB_TE') + 10.^(0.1*acoustic.LpA_Total));
                end
                if (method.BE==1) 
                for i = 1 : geom.Nb
                   acoustic.BE.OASPLA_Total(i) = 10*log10(10^(0.1*acoustic.BE.OASPLA_Tonal(i))+10^(0.1*acoustic.BE.OASPLA_BB_TE(i)));
                end
                end
            elseif (method.broadband==2)
                if (method.optim==-1)
                    fprintf('%s\n','Broadband noise computation (turbulence ingestion)');
                end
                acoustic            = Acoustic_broadband_TI(test,geom,flow,acoustic,obs,data,method); %broadband noise through Rozenberg et al's method
                %
                if (method.freqband==0)
                   acoustic.OASPL_tot    = 10*log10(10^(0.1*acoustic.OASPL_Total)+10^(0.1*acoustic.OASPL_BB_TI));
                   acoustic.Lp_Total     = 10*log10(10.^(0.1*acoustic.Spp_BB_TI') + 10.^(0.1*acoustic.Lp_Total));
                elseif (method.freqband==1)
                   acoustic.OASPLA_Total = 10*log10(10^(0.1*acoustic.OASPLA_Total)+10^(0.1*acoustic.OASPLA_BB_TI));
                   acoustic.LpA_Total    = 10*log10(10.^(0.1*acoustic.SppA_BB_TI') + 10.^(0.1*acoustic.LpA_Total));
                end
                if (method.BE==1) 
                for i = 1 : geom.Nb
                   acoustic.BE.OASPLA_Total(i) = 10*log10(10^(0.1*acoustic.BE.OASPLA_Tonal(i))+10^(0.1*acoustic.BE.OASPLA_BB_TI(i)));
                end
                end
            elseif (method.broadband==3)
                if (method.optim==-1)
                    fprintf('%s\n','Broadband noise computation (vortex shedding)');
                end
            elseif (method.broadband==4)
                if (method.optim==-1)
                    fprintf('%s\n','Broadband noise computation (all sources)');
                end
            elseif (method.broadband==5)
                if (method.optim==-1)
                    fprintf('%s\n','Broadband noise computation (separation)');
                end
                [flow,acoustic]         = Acoustic_separation(test,geom,flow,acoustic,obs,data,method);
                
                acoustic.OASPL_Total    = 10*log10(10^(0.1*acoustic.OASPL_Total)+10^(0.1*acoustic.OASPL_sep));
                Spp_sep                 = interp1(acoustic.fbb,acoustic.Spp_sep,acoustic.f);
                k                       = 1;
                while isnan(Spp_sep(k)) == 1 % avoid NaN from interpolation
                    k                   = k+1;
                end
                if k > 1
                    Spp_sep(1:k-1)      = Spp_sep(k);
                end
                k                       = length(Spp_sep);
                while isnan(Spp_sep(k)) == 1
                    k                   = k-1;
                end
                if k < length(Spp_sep)
                    Spp_sep(k+1:end)    = Spp_sep(k);
                end
                acoustic.Spp_Total      = 10*log10(10.^(0.1*Spp_sep) + 10.^(0.1*acoustic.Spp_Total));
                if method.dBA == 1
                    acoustic.OASPLA_Total = 10*log10(10^(0.1*acoustic.OASPLA_Total)+10^(0.1*acoustic.OASPLA_sep));
                    Spp_A_sep           = interp1(acoustic.fbb,acoustic.Spp_A_sep,acoustic.f);
                    k                   = 1;
                    while isnan(Spp_A_sep(k)) == 1 % avoid NaN from interpolation
                        k = k+1;
                    end
                    if k > 1
                        Spp_A_sep(1:k-1) = Spp_A_sep(k);
                    end
                    k                   = length(Spp_A_sep);
                    while isnan(Spp_A_sep(k)) == 1
                        k               = k-1;
                    end
                    if k < length(Spp_A_sep)
                        Spp_A_sep(k+1:end) = Spp_A_sep(k);
                    end
                    acoustic.Spp_A_Total  = 10*log10(10.^(0.1*Spp_A_sep) + 10.^(0.1*acoustic.Spp_A_Total));
                end
            end
        else
            fprintf('%s\n','Acoustic module not called : non convergence of iso-thrust');
        end
    else
        fprintf('%s\n','Acoustic module not called : blade too stalled (eliminated)');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
