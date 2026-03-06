function [Polar] = Aerodynamic_polarfunction(polar_file)

%---- read f( Cl , Cd ) function file

% input :
% polar_file     Name of the Polar file

% output :
% Polar data structure
% Polar.Name      = 'optim';
% Polar.Airfoil   = 'goe265';
% Polar.ReRef     = 10 000;
% Polar.ReExp     = 10 000;
    
fid = fopen( polar_file , 'r' ) ;
LiftFunction = textscan(fid, '%f %f %f %f', 'HeaderLines', 1) ;
Polar.CL0     = LiftFunction{ 1 } ;
Polar.DCLDA   = LiftFunction{ 2 } ;
Polar.CLMIN   = LiftFunction{ 3 } ;
Polar.CLMAX   = LiftFunction{ 4 } ;
fclose( fid ) ;

fid = fopen( polar_file , 'r' ) ;
DragFunction = textscan(fid, '%f %f %f %f', 'HeaderLines', 4) ;
Polar.CD0     = DragFunction{ 1 } ;
Polar.CD2U    = DragFunction{ 2 } ;
Polar.CD2L    = DragFunction{ 3 } ;
Polar.CLCD0   = DragFunction{ 4 } ;
fclose( fid ) ;

fid = fopen( polar_file , 'r' ) ;
ReynoldsFunction = textscan(fid, '%f %f %f %f', 'HeaderLines', 7) ;
Polar.REREF = ReynoldsFunction{ 1 } ;
Polar.REEXP = ReynoldsFunction{ 2 } ;
fclose( fid ) ;
