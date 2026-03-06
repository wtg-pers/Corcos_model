%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary Layer informations Loader
%
% WARNING! values are non-dimensionnal !
%
% V1.0 by R.SERRE and N.PICHON 23/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [boundary] = Aerodynamic_MonoAirfoil_BoundaryLayer(data,geom)

% AIRFOIL SECTION ===========================================================
FOIL = load(['./AirfoilDataBase/',geom.airfoil,'.dat']);
Xp=FOIL(:,1);Np=length(Xp);
Yp=FOIL(:,2);
% BLADE THICKNESS ===========================================================
% Position pour blade thickess
id_ext=20;
id_int=Np-id_ext;
%
boundary.Dt = abs(Yp(id_ext)-Yp(id_int));


% load data '.mat' file
load(['./BoundaryLayerDataBase/',geom.airfoil,'.mat']);

% read data
for i = 1 : size(BLayer,2) % for each Reynolds number
    for j = 1 : length(BLayer(i).Alpha) % for each Angle of Attack
        %
        % XFOIL boundary informations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rc          = BLayer(i).Re;
        AOA         = BLayer(i).Alpha(j);
        Vinf        = Rc*data.nu ;
        Ns          = numel(BLayer(i).BoundaryLayer(j).s);
        x           = zeros(Ns,1);
        y           = x;
        Dstar       = x;
        D           = x;
        Theta       = x;
        UeVinf      = x;
        Cf          = x;
        %
        for k = 1 : Ns %for each point of the airfoil
            s               = BLayer(i).BoundaryLayer(j).s(k);
            x(k)            = BLayer(i).BoundaryLayer(j).x(k);
            y(k)            = BLayer(i).BoundaryLayer(j).y(k);
            UeVinf(k)       = BLayer(i).BoundaryLayer(j).UeVinf(k);
            Dstar(k)        = BLayer(i).BoundaryLayer(j).Dstar(k);
            Theta(k)        = BLayer(i).BoundaryLayer(j).Theta(k);
            Cf(k)           = BLayer(i).BoundaryLayer(j).Cf(k);
            H               = BLayer(i).BoundaryLayer(j).H(k);
            % BL-edge mach number
            Me              = UeVinf(k)  * Vinf / data.c ;
            % Drela, Giles, AIAAj VOL.25 NO.10, equation (23)
            Hk              = ( H - 0.29*Me^2 )/( 1 + 0.113*Me^2 );
            D(k)            = Theta(k)*(3.15 + (1.72/(Hk-1)))+Dstar(k);
            
            if ( (x(k) == x(1)) && (y(k) == -y(1)) )
                TEp         = k-4;
                TEs         = 1+4;
            end
            
        end
        
        %export data
        boundary.Uep(j,i)  = -UeVinf(TEp)*Vinf;
        boundary.Ues(j,i)  = UeVinf(TEs)*Vinf;
        boundary.Cfp(j,i)  = Cf(TEp);
        boundary.Cfs(j,i)  = Cf(TEs);
        boundary.tp(j,i)   = Theta(TEp);
        boundary.ts(j,i)   = Theta(TEs);
        boundary.dsp(j,i)  = Dstar(TEp);
        boundary.dss(j,i)  = Dstar(TEs);
        boundary.dp(j,i)   = D(TEp);
        boundary.ds(j,i)   = D(TEs);
        boundary.Re(i)     = Rc;
        boundary.Alpha(j)  = AOA;
    end 
end


for i = 1 : size(boundary.Uep,2);
    j0       = size(boundary.Uep,1);
    while  boundary.Uep(j0,i) == 0
        j0   = j0-1;
    end
    
    for j = min(j0 + 1,size(boundary.Uep,1)) : size(boundary.Uep,1)
        boundary.Uep(j,i)  = boundary.Uep(j0,i);
        boundary.Ues(j,i)  = boundary.Ues(j0,i);
        boundary.Cfp(j,i)  = boundary.Cfp(j0,i);
        boundary.Cfs(j,i)  = boundary.Cfs(j0,i);
        boundary.tp(j,i)   = boundary.tp(j0,i);
        boundary.ts(j,i)   = boundary.ts(j0,i);
        boundary.dsp(j,i)  = boundary.dsp(j0,i);
        boundary.dss(j,i)  = boundary.dss(j0,i);
        boundary.dp(j,i)   = boundary.dp(j0,i);
        boundary.ds(j,i)   = boundary.ds(j0,i);
        boundary.Alpha(j)  = boundary.Alpha(j0);
    end
end
