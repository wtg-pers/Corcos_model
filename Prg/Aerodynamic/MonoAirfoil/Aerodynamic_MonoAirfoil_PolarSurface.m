function [polar] = Aerodynamic_MonoAirfoil_PolarSurface(geom,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polar surface creation
%
% V1.0 by N.PICHON 03/10/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%look up available reynolds for concerned airfoil
airfoil = geom.airfoil;
if (method.fitpolar==1)
polarpath = ['./PolarDataBase/' airfoil 'Fitted'];
else
polarpath = ['./PolarDataBase/' airfoil];
end
Re      = load([polarpath '/Reynolds.dat']);
%
%create AoA vector
Alpha   = load([polarpath '/Alpha.dat']);

%initialise Lift and Drag coefficients
Cl      = zeros(length(Re),length(Alpha))';
Cd      = zeros(length(Re),length(Alpha))';

for i = 1 : length(Re) %for each Reynolds number
    aa      = load([polarpath '/Re',num2str(Re(i)),'_Alpha.dat']);  %load up the available alphas vector
    if isequal(aa,Alpha)==1 %if this vector is complete (i.e. = Alpha), just load Cl and Cd tables        
        a   = load([polarpath '/Re',num2str(Re(i)),'_Cl.dat']);
        Cl(1:length(aa),i)  = a(1:length(aa));
        a   = load([polarpath '/Re',num2str(Re(i)),'_Cd.dat']);
        Cd(1:length(aa),i)  = a(1:length(aa));
        
    else % otherwise interpolate missing values
        a           = load([polarpath '/Re',num2str(Re(i)),'_Cl.dat']);
        b           = a(1:length(aa));
        Cl(:,i)     = interp1(aa,b,Alpha,'spline');
        
        a           = load([polarpath '/Re',num2str(Re(i)),'_Cd.dat']);
        b           = a(1:length(aa));
        Cd(:,i)     = interp1(aa,b,Alpha,'spline');
        
        if aa(end)<45 % forbid extrapolation. A cap is put instead
            [~,j]           = min(abs(aa(end)-Alpha));
            Cl(j+1:end,i)   = Cl(j,i) * ones(1,length(Cl)-j);
            Cd(j+1:end,i)   = Cd(j,i) * ones(1,length(Cd)-j);
        end        
    end
end

% record data
polar.Cl        = Cl;
polar.Cd        = Cd;
polar.Alpha     = Alpha;
polar.Re        = Re;
