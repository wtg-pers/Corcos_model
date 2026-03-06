function [flow,geom] = Aerodynamic_MultiAirfoil_BEMT(test,geom,polar,data,method)

Nblades         = geom.NBlades;

BladePitchAngle = 0;    % Equivalent to the DBE variable;

%refine spatial ( radial ) resolution
[R,C,B,Fx,N]    = Geom_boundary (geom);
RAD             = R(end);
Rhub            = R(1);

%% general parameters (harcoded values)
VSO     = data.c ;          % 'VS0' is the speed of sound ( m / s )
RHO     = data.rho ;        % 'RHO' is the fluid density ( kg / m3 )
RMU     = data.mu ;      % 'RMU' is the fluid dynamic viscosity ( kg / m / s )


%% Actual BEMT resolution phase
VEL     = test.Speed;
OMG     = test.RPM * 2 *pi / 60 ;                                      % 'OMG' is the rotationnal velocity ( rad / s )
DBE     = BladePitchAngle * pi / 180 ;                                    % 'DBE' is the collective pitch angle ( rad )

%% go over radial stations

% Pre allocation
dTdr    = zeros(1,N);
dQdr    = zeros(1,N);
A       = dTdr;
MA      = dTdr;
Vr      = dTdr;
a       = dTdr;
b       = dTdr;
Vt      = dTdr;

%SA_AOA = dlmread('./Prg/Aerodynamic/MultiAirfoil/AngleOfAttackSingleAirfoil.dat');SA_AOA=SA_AOA*pi/180;
%localAlpha = 0;

for j = 1 : N
    bladeSolidityFactor = Nblades*C(j)/(2*pi*R(j));
    localGeometricSpeed = sqrt(VEL*VEL + (OMG*R(j))^2);
    
    localReynolds = localGeometricSpeed*C(j)*RHO/RMU;
    
    
    
    % First guesses for the Bisection method
    negativeInflowAngle = 0;
    positiveInflowAngle = pi/2;
    
%while ( abs(localAlpha - SA_AOA(j))*180/pi > 0.01)
%fprintf(' Ecart : %2.2f \n',abs(localAlpha - SA_AOA(j))*180/pi)
    % Bisection method to solve the transcendental equation
    for n = 1 : 15
        
        newInflowAngle = 0.5*(negativeInflowAngle + positiveInflowAngle);
        localAlpha = DBE + B(j) - newInflowAngle;
        
        
%%%%%%%%%%%%%calculation of local lift and drag coefficients%%%%%%%%%%%%%%%
        if method.polar == 1
            if (method.NFoilOpt==N)
            airfoilindex=j;
            else
            airfoilindex=ceil(j/(N/method.NFoilOpt));
            end
            [Cd_loc,Cl_loc] = Aerodynamic_MultiAirfoil_PolarInterpolation(localAlpha,localReynolds,polar,airfoilindex);
        else
            [Cl_loc,Cd_loc] = Aerodynamic_CLCDFUN(polar,localReynolds,localAlpha,VSO,RMU,RHO);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tanGama = Cd_loc/Cl_loc;
        gama = atan(tanGama);
        
        % Non stationary case
        if VEL~=0
            secGama = sec(gama);
            newFunctionValue = (OMG*R(j)*sin(newInflowAngle) ...
                - VEL*cos(newInflowAngle))*sin(newInflowAngle) ...
                - 0.25*bladeSolidityFactor*Cl_loc*secGama*(OMG*R(j)*cos(newInflowAngle+gama) ...
                + VEL*sin(newInflowAngle + gama));
            
            if(newFunctionValue < 0)
                negativeInflowAngle = newInflowAngle;
            else
                positiveInflowAngle = newInflowAngle;
            end
            % Stationary case
        else
            newFunctionValue = 4*(sin(newInflowAngle))^2 - bladeSolidityFactor*Cl_loc*sqrt(1+tanGama^2)*cos(newInflowAngle + gama);
            if(newFunctionValue < 0)
                negativeInflowAngle = newInflowAngle;
            else
                positiveInflowAngle = newInflowAngle;
            end
            
        end
           
    end     
%% Adjusting Twist angle to fit with angle of attack of equivalent single airfoil version
%B(j) = B(j) - (localAlpha-SA_AOA(j));
%end
        % Non stationary case
        if VEL~=0
            F = (sin(newInflowAngle) - 0.25*bladeSolidityFactor*Cl_loc*sec(gama)*csc(newInflowAngle)*cos(newInflowAngle + gama));
            G = OMG*R(j)*F/VEL;
            Vr(j) = VEL/F;
            a(j) = sin(newInflowAngle)/F -1;
            b(j) = 1 - cos(newInflowAngle)/G;
        else % Stationary case
            if (Cl_loc==0)
                b(j) = 0;
            else
                b(j) = 1/(1 + 4*cos(gama)*sin(newInflowAngle)*cos(newInflowAngle)/(bladeSolidityFactor*Cl_loc*sin(newInflowAngle + gama)));
            end
            Vr(j) = (1-b(j))*OMG*R(j)*sec(newInflowAngle);
            Vt(j) = Vr(j)*sin(newInflowAngle);
        end
%%%%%%%%%%%%%%%%Prandtl's tip and hub loss model %%%%%%%%%%%%%%%%%%%%%%%%%%
        if method.tip == 0
            Ftip = 1.0;
        elseif method.tip == 1
            f = (Nblades/2)*(RAD/R(j) - 1)*(1/sin(B(end)));
            if B(end) == 0 && RAD == R(j)
                f = 0;
            end
            Ftip = (2/pi) * acos(exp(-f));
        end
        if method.hub == 0
            Fhub = 1.0;
        elseif method.hub == 1
            f = (Nblades/2)*(1-Rhub/R(j))*(1/sin(B(1))); 
            if B(1) == 0 && Rhub == R(j)
                f = 0;
            end
            Fhub = (2/pi) * acos(exp(-f));
        elseif method.hub == 2
            f = (Nblades/2)*(1-Rhub/R(j))*(1/sin(B(1)));
            f_RAD=(Nblades/2)*(1-Rhub/RAD)*(1/sin(B(1)));
            if B(1) == 0 && Rhub == R(j)
                f = 0;
            end
            Fhub = acos(exp(-f))/acos(exp(-f_RAD));
        elseif method.hub == 3
            f = (Nblades/2)*(R(j)/Rhub-1)*(1/sin(B(1)));
            if B(1) == 0 && Rhub == R(j)
                f = 0;
            end            
            Fhub = (2/pi) * acos(exp(-f));
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Elemental thrust and torque:
        if VEL~=0
            dTdr(j) = Ftip*Fhub*4*pi*RHO*R(j)*VEL*VEL*a(j)*(1+a(j));
            dQdr(j) = Ftip*Fhub*4*pi*RHO*(R(j)^3)*OMG*VEL*(1+a(j))*b(j);
        else
            dTdr(j) = Ftip*Fhub*4*pi*RHO*R(j)*((Vr(j)*sin(newInflowAngle))^2);
            dQdr(j) = Ftip*Fhub*4*pi*RHO*(R(j)^3)*b(j)*OMG*Vr(j)*sin(newInflowAngle);
        end
        
        MA( j ) = OMG*R(j)/VSO;
        A( j )  = localAlpha;
        
%%%%%%%%% sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if method.sweep == 1
            delta=atan(Fx(j)/R(j));
            
            dTdr(j)=dTdr(j)*cos(delta);
            dQdr(j)=dQdr(j)*cos(delta)^2;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%outputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
flow.newFunctionValue(j) = newFunctionValue;
flow.AOA(j)             = A(j)*180/pi;
flow.Lift(j)            = dTdr(j)         / Nblades; % N/m (normalized by Blade ELement)
flow.Drag(j)            = dQdr(j) /  R(j) / Nblades; 
flow.LocalReynolds(j)   = localReynolds;
flow.Cl_local(j)        = Cl_loc;
flow.Cd_local(j)        = Cd_loc;
flow.LocalThrust(j)     = dTdr(j)        / Nblades;
flow.LocalTorque(j)     = dQdr(j) / R(j) / Nblades;
flow.Mach(j)            = MA(j);
flow.InducedVelocity(j) = Vr(j); % La vitesse induite est de 10% de la vitesse de rotation (Thierry Jardin)
flow.TangentialSpeed(j) = Vt(j);

end

flow.Thrust             = trapz(R,dTdr);
flow.Torque             = trapz(R,dQdr);
geom.Rb                 = R;
geom.Cb                 = C;
geom.Tb                 = B*180/pi;
geom.Fxb                = Fx;
geom.Nb                 = N;
%
%geom.T = spline(geom.Rb,geom.Tb,geom.Rq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geom_Inverseboundary(geom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%figure
%plot(flow.InducedVelocity)
%title('NP')
end
