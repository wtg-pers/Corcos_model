function [polar] = Aerodynamic_MultiAirfoil_PolarSurface(geom,method)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polar surface creation
%
% V1.0 by N.PICHON 03/10/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /!\ :
%look up available reynolds for concerned airfoil
%
Re      = dlmread(['./AirfoilOptimization/PolarDataBase/Reynolds.dat']);
Alpha   = dlmread(['./AirfoilOptimization/PolarDataBase/Alpha.dat']);
%
for k = 1 : method.NFoilOpt
   wd = ['./AirfoilOptimization/PolarDataBase/Airfoil' sprintf('%03d',k) '/'];
   %
   %initialise Lift and Drag coefficients
   Cl      = zeros(length(Re),length(Alpha))';
   Cd      = zeros(length(Re),length(Alpha))';
   %
   for i = 1 : length(Re) %for each Reynolds number
       aa      = load([wd 'Re',num2str(Re(i)),'_Alpha.dat']);  %load up the available alphas vector
       if isequal(aa,Alpha)==1 %if this vector is complete (i.e. = Alpha), just load Cl and Cd tables        
           a   = load([wd 'Re',num2str(Re(i)),'_Cl.dat']);
           Cl(1:length(aa),i)  = a(1:length(aa));
           a   = load([wd 'Re',num2str(Re(i)),'_Cd.dat']);
           Cd(1:length(aa),i)  = a(1:length(aa));
           
       else % otherwise interpolate missing values
           a           = load([wd 'Re',num2str(Re(i)),'_Cl.dat']);
           b           = a(1:length(aa));
           Cl(:,i)     = interp1(aa,b,Alpha,'spline');
           
           a           = load([wd 'Re',num2str(Re(i)),'_Cd.dat']);
           b           = a(1:length(aa));
           Cd(:,i)     = interp1(aa,b,Alpha,'spline');
           
           if aa(end)<45 % forbid extrapolation. A cap is put instead
               [~,j]           = min(abs(aa(end)-Alpha));
               Cl(j+1:end,i)   = Cl(j,i) * ones(1,length(Cl)-j);
               Cd(j+1:end,i)   = Cd(j,i) * ones(1,length(Cd)-j);
           end        
       end
   end
   %
   % record data
   polar.airfopt(k).Cl        = Cl;
   polar.airfopt(k).Cd        = Cd;
   polar.airfopt(k).Alpha     = Alpha;
   polar.airfopt(k).Re        = Re;
end
