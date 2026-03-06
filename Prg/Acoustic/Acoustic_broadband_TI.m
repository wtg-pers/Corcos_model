%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbulence Interaction Noise
% Model from Roger Moreau IJA 2010
%
% V1.0 by R.SERRE
% 16/05/2017 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [acoustic] = Acoustic_broadband_TI(test,geom,flow,acoustic,obs,data,method)
%
% Load Boundary layer informations
if (method.airfopt==0)
boundary    = Aerodynamic_MonoAirfoil_BoundaryLayer(data,geom);
end
%
%Initialisation
Spp         = zeros(length(acoustic.fbb),1);
Spp_k       = zeros(length(geom.Rb),length(acoustic.fbb));
Doppler     = zeros(acoustic.Ntobs,1);
Spp_psi     = zeros(acoustic.Ntobs,1);
PSI         = [0:acoustic.dtobs:acoustic.Tobs]*test.RPM*2*pi/60;    %angular step
%
sinTHETA    = obs.x0(1)/obs.Robs;           % theta is the azimuthal angle
%
%calculations
%%for k = round( (length(geom.Rb)+1)/2) :  length(geom.Rb)            %for each blade element
%for k = round( (3*(length(geom.Rb)+1))/4) :  round( (3*(length(geom.Rb)+1))/4)            %for each blade element
%for k = length(geom.Rb) : length(geom.Rb)             %for each blade element
%for k = round( (3*(length(geom.Rb)+1))/4) : length(geom.Rb)             %for each blade element
for k = 1 : length(geom.Rb)             %for each blade element
    %
    if (method.airfopt==1)
       if (method.NFoilOpt==length(geom.Rb))
          airfoilindex=k;
       else
          airfoilindex=ceil(k/(length(geom.Rb)/method.NFoilOpt));
       end
       boundary    = Aerodynamic_MultiAirfoil_BoundaryLayer(data,geom,method,airfoilindex);
    end
    %
    %init values relative to blade element
    Radius   = geom.Rb(k);                   % Rayon de l'element de pale
    Chord    = geom.Cb(k);                   % Chord de l'element de pale
    U        = (test.RPM*2*pi/60) * Radius;  % blade element speed
    Rc       = U*Chord / data.nu;            % local Reynolds
    Mach     = U/data.c;                     % Mach number
    Mt       = Mach;                         % Rotational Mach number
    b2       = 1 - (Mach)^2;                 % Prandtl-Glauert Compressibility coefficient
    L        = geom.dR;                      % Blade element span-wise length
    %
    % load corresponding boundary layer values
    [~,indexOfClosestReynolds] = min(abs(boundary.Re-Rc)); % find closest available index in boundary layer tables
    [~,indexOfLocalAlpha]      = min(abs(boundary.Alpha - flow.AOA(k)));
    %
    dp  = boundary.dp(indexOfLocalAlpha,indexOfClosestReynolds)  * Chord ;
    ds  = boundary.ds(indexOfLocalAlpha,indexOfClosestReynolds)  * Chord ;
    dsp = boundary.dsp(indexOfLocalAlpha,indexOfClosestReynolds) * Chord ;
    dss = boundary.dss(indexOfLocalAlpha,indexOfClosestReynolds) * Chord ;
    Uep = boundary.Uep(indexOfLocalAlpha,indexOfClosestReynolds);% / Chord ; % Modified 16/05/2017
    Ues = boundary.Ues(indexOfLocalAlpha,indexOfClosestReynolds);% / Chord ;
    Cfp = boundary.Cfp(indexOfLocalAlpha,indexOfClosestReynolds);
    Cfs = boundary.Cfs(indexOfLocalAlpha,indexOfClosestReynolds);
    
    %
    if (method.freqband==0)
       startf=2;
    elseif (method.freqband==1)
       startf=1;
    end
    %
    % Integration in transverse wave number 
   % NWN   = 11;
   % WNmin = -10;
    %WNmax =  10;
    %dWN   = (WNmax-WNmin)/(NWN-1);
   % WNV   = [WNmin : dWN : WNmax];
    %
    for i=startf:length(acoustic.fbb) % for each frequency band
        wobs = 2*pi*acoustic.fbb(i); %associated pulsation (from observer point of view
        for j=1:acoustic.Ntobs %for each timestep
            Doppler(j) = 1 + Mt*sin(PSI(j))*sinTHETA - Mt*cos(PSI(j))*obs.x0(2)/obs.Robs; % Doppler factor due to blade rotation
            % observer's position in the moving frame
            Y1 = obs.x0(1) - geom.Rb(k)*cos(PSI(j));
            Y2 = obs.x0(2) - geom.Rb(k)*sin(PSI(j));
            Y3 = obs.x0(3);
            S0 = sqrt( Y1^2 + b2*(Y2^2 + Y3^2) );
            %
            wsrc = Doppler(j)*wobs; % Source actual emission pulsation
            %
            Fintg=0;
            kac    =  wsrc / data.c;   
            k1     =  wsrc /  U ;
            ks1    = k1 * (Chord/2) ;
            
            %%%%%%%%%%%%%%%% Modified by Viswesh SUJJUR BALARAMRAJA 13/6/2017%%%%%%%%%%%%
            ks2min =  -(2*L/Chord)*kac*Chord/(2*sqrt(b2));
            ks2max =  (2*L/Chord)*kac*Chord/(2*sqrt(b2));
            dks2   = 0.05*kac*Chord/(2*sqrt(b2));
            dk2   =  dks2/(Chord/2);
           
   
            for ks2 = ks2min:dks2: ks2max ;
            % LTI ===============================
            %
            kac    =  wsrc / data.c;
            %k1     =  wsrc /  U ;
            k2     =  ks2 / (Chord/2);
            mub    = (kac*Chord) / (2*b2) ;
            %ks1    = k1 * (Chord/2) ;
            %ks2    = k2 * (Chord/2);
            kappa  = sqrt( (ks2^2/b2) - mub^2 );
            T0     = ks1 * Mach / (sqrt(b2)*ks2);
            T2     = mub *(Mach - Y1*S0) - 0.25*pi;
            T3     = kappa + 1i*mub*(Y1/S0) ;
            T4     = 1i*kappa - mub*(Y1/S0);
            kappa2 = mub^2 * ((1/T0^2)-1.0);
            %
            % TURBULENCE SPECTRUM - Amiet JSV 1975 
            %
            TurbLength = boundary.Dt*Chord + dsp + dss ; % Wake Width matching Taylor Scale of Turbulence 
            %TurbLength = boundary.Dt*Chord + dp + ds ; % Wake Width matching Integrale Scale of Turbulence 
            %VelScale   = 0.5*(Uep + Ues) ;    % Order of magnitude of longitudinal velocity
            VelScale = mean(flow.TangentialSpeed(:));
            %
            ke = sqrt(pi)/TurbLength * gamma(5./6.) / gamma(1./3.) ;
            %
            kh1     = k1 / ke ;
            kh2     = k2 / ke ;
            
            % von Karman model : 
            PHIww = 4/(9*pi) * (VelScale/ke)^2 * (kh1^2 + kh2^2)/(1 + kh1^2 + kh2^2)^(7/3) ;
            % Liepmann model : 
            %PHIww = (3*VelScale^2)/(4*pi) * TurbLength^2 * (TurbLength^2 * (k1^2 + k2^2))/( 1 + TurbLength^2 * (k1^2 + k2^2))^(5/2);
            % 
            % Fresnel Integrals Discretisation
            nt=101;
            % SUB-CRITICAL GUSTS ----------------
            if (T0 < 1)
               %
               %xi    = (2*1i*T3);tt=[0:xi/(nt-1):xi];
               %%xi    = real(2*1i*T3);tt=[0:xi/(nt-1):xi];
               %Fresnel1 = 0;
               %for jj=2:nt-1
               %   Fresnel1 = Fresnel1 + exp( 1i*tt(jj) ) / sqrt( 2*pi*tt(jj) ) * 0.5*( tt(jj+1) - tt(jj-1) );
               %end
   [cfren1,sfren1] = fresnelMoser(sqrt(2*1i*T3));Fresnel1=cfren1+1i*sfren1;
   %[cfren1,sfren1] = fresnelMoser(2*1i*T3);Fresnel1=sqrt(cfren1^2+sfren1^2);
               xi    = sqrt(4*abs(kappa));tt=[0:xi/(nt-1):xi];
               % fprintf('%f\n', xi);
               % stop
               % for aa=1:1:nt
               %    fprintf('%i %f\n', aa,tt(aa)); 
               % end
               
               % stop
               ferr = 0;
               for jj=2:nt-1
                  ferr = ferr + exp( -tt(jj)^2 ) * 0.5*( tt(jj+1) - tt(jj-1) );
               end
               ferr=2/sqrt(pi) * ferr;
               %xi    = (2*( 1i*kappa + mub*(Y1/S0)));tt=[0:xi/(nt-1):xi];
               %%xi    = real(2*( 1i*kappa + mub*(Y1/S0)));tt=[0:xi/(nt-1):xi];
               %Fresnel2 = 0;
               %for jj=2:nt-1
               %   Fresnel2 = Fresnel2 + exp( 1i*tt(jj) ) / sqrt( 2*pi*tt(jj) ) * 0.5*( tt(jj+1) - tt(jj-1) );
               %end
   [cfren2,sfren2] = fresnelMoser(sqrt(2*( 1i*kappa + mub*(Y1/S0))));Fresnel2=cfren2+1i*sfren2;
   %[cfren2,sfren2] = fresnelMoser((2*( 1i*kappa + mub*(Y1/S0))));Fresnel2=sqrt(cfren2^2+sfren2^2);
               %
               L1 = (-1/pi)*sqrt( 2 / ( (ks1 + 1i*b2*kappa)*1i*T3 ) )*exp(-1i*T2)*Fresnel1;
               %
               L2 = ( (exp(-1i*T2))/(pi*sqrt(2*pi*(ks1+1i*b2*kappa))*T3) ) ...
                  * ( 1 - exp(-2*T3) - ferr + 2*exp(-2*T3)*sqrt(kappa/(1i*kappa + mub*(Y1/S0)))*Fresnel2 );
               %
            % SUPERCRITICAL GUSTS ---------------
            elseif (T0 > 1 && kappa2 < 0)
               %
               %xi    = 2*T4;tt=[0:xi/(nt-1):xi];
               %Fresnel3 = 0;
               %for jj=2:nt-1
               %   Fresnel3 = Fresnel3 + exp( 1i*tt(jj) ) / sqrt( 2*pi*tt(jj) ) * 0.5*( tt(jj+1) - tt(jj-1) );
               %end
   [cfren3,sfren3] = fresnelMoser(sqrt(2*T4));Fresnel3=cfren3+1i*sfren3;
   %[cfren3,sfren3] = fresnelMoser(2*T4);Fresnel3=sqrt(cfren3^2+sfren3^2);
               %xi    = 4*1i*kappa;tt=[0:xi/(nt-1):xi];
               %Fresnel4 = 0;
               %for jj=2:nt-1
               %   Fresnel4 = Fresnel4 + exp( 1i*tt(jj) ) / sqrt( 2*pi*tt(jj) ) * 0.5*( tt(jj+1) - tt(jj-1) );
               %end
   [cfren4,sfren4] = fresnelMoser(sqrt(4*1i*kappa));Fresnel4=cfren4+1i*sfren4;
   %[cfren4,sfren4] = fresnelMoser(4*1i*kappa);Fresnel4=sqrt(cfren4^2+sfren4^2);
               %xi    = 2*(1i*kappa + mub*(Y1/S0));tt=[0:xi/(nt-1):xi];
               %Fresnel5 = 0;
               %for jj=2:nt-1
               %   Fresnel5 = Fresnel5 + exp( 1i*tt(jj) ) / sqrt( 2*pi*tt(jj) ) * 0.5*( tt(jj+1) - tt(jj-1) );
               %end
   [cfren5,sfren5] = fresnelMoser(sqrt(2*( 1i*kappa + mub*(Y1/S0))));Fresnel5=cfren5+1i*sfren5;
   %[cfren5,sfren5] = fresnelMoser((2*( 1i*kappa + mub*(Y1/S0))));Fresnel5=sqrt(cfren5^2+sfren5^2);
               %
               L1 = (-1/pi)*sqrt( 2 / ( (ks1 + 1i*b2*kappa)*1i*T4 ) )*exp(-1i*T2)*Fresnel3;
               %
               L2 = ( (exp(-1i*T2))/(pi*sqrt(2*pi*(ks1+b2*1i*kappa))*T4) ) ...
                  * ( 1i*(1 - exp(2*1i*T4)) - (1-1i)*( Fresnel4 - exp(2*1i*T4)*sqrt(2*1i*kappa/(1i*kappa + mub*(Y1/S0)))*Fresnel5 ));
               %
            else
               fprintf(' GUST PROBLEM \n')
               L1 = 0; L2 = 0;
            end
            
            % ===================================
            LTI = L1 + L2 ;
            % ===================================
            %
            Fcorr = (sin( (kac*Y2/S0 - k2)*(L/2) )^2)/( pi * (L/2) * (kac*Y2/S0 - k2)^2 ) ; 
            %
            Fintg = Fintg + (PHIww * abs(LTI)^2 * Fcorr)*dk2 ;
            %
            end 
         %
         Spp_psi(j) = ((data.rho*kac*Chord*Y3)/(2*S0^2))^2 * pi*U*(L/2) * Fintg;
         %
        end
        Spp_k(k,i)   = Spp_k(k,i) + geom.NBlades * trapz(Doppler.*Spp_psi)/(2*pi); %time integration
        
        Spp(i) = (sqrt(Spp(i)) + sqrt(abs(Spp_k(k,i))))^2; % sum over blade elements
    end
end
%
%Overall trailing edge noise calculation
OASPL_BB_TI = 0;
if (method.freqband==0)
    bw = ones (1,length(acoustic.fbb)); % bandwith is constant equals to 1Hz if tonal noise frequency spacing was selected
else
    bw = zeros(1,length(acoustic.fbb));
    for i = 1 : length(acoustic.fbb)-1 % determines bandwidth when freq bands are equally spaced in the log space
        bw(i) = 2*(10.^((log10(acoustic.fbb(i+1))+log10(acoustic.fbb(i)))/2)-acoustic.fbb(i));
    end
    bw(length(acoustic.fbb)) = 2*(-10.^((log10(acoustic.fbb(end))+log10(acoustic.fbb(end-1)))/2)+acoustic.fbb(end));
end

%
Spp   = Spp.*bw';
for i = 1 : geom.Nb
Spp_k(i,:) = Spp_k(i,:).*bw(:)';
end
%
for i = 1 : length(acoustic.fbb) %sum for all frequency bands (multiplied by bandwith)
    OASPL_BB_TI = OASPL_BB_TI + Spp(i);
end
%OASPL_BB_TI
acoustic.OASPL_BB_TI = 10*log10(OASPL_BB_TI/(data.pref^2));
%
Spp_A=zeros(length(acoustic.fupper),1);
%A-weighting
if (method.freqband==0)
   k=1;
   for j=1:length(acoustic.fbb)
      if (round(acoustic.fbb(j)) > acoustic.flower(k) && round(acoustic.fbb(j)) < acoustic.fupper(k))
        Spp_A(k)  = Spp_A(k) + (data.pref^2)*10.^( log10(Spp(j)/(data.pref^2)) );
      elseif (round(acoustic.fbb(j))==acoustic.fupper(k))
        Spp_A(k)  = 10.0*log10(Spp_A(k)/data.pref^2) + acoustic.ponderation(k);
         k=k+1;
      end
   end
    OASPLA_BB_TI = 0;
    for i = 1 : length(Spp_A) %sum for all frequency bands (multiplied by bandwith)
        OASPLA_BB_TI = OASPLA_BB_TI + Spp_A(i);
    end
    acoustic.OASPLA_BB_TI   = 10*log10(OASPLA_BB_TI/(data.pref^2));
elseif (method.freqband==1)
    OASPLA_BB_TI = 0;
    for i = 1 : length(acoustic.fbb) %sum for all frequency bands (multiplied by bandwith)
        Spp_A(i)   = (data.pref^2)*10.^( log10(Spp(i)/(data.pref^2))+acoustic.ponderation(i)/10 );
        OASPLA_BB_TI = OASPLA_BB_TI + Spp_A(i);
    end
    acoustic.OASPLA_BB_TI   = 10*log10(OASPLA_BB_TI/(data.pref^2));
end
%export data
acoustic.Spp_BB_TI  = 10*log10(Spp/(data.pref^2));
acoustic.SppA_BB_TI = 10*log10(Spp_A/(data.pref^2));
%
%Blade element post-processing
%% CORRECT WITH BANDWITH AND A-WEIGHTING
%if method.BE == 1
%    acoustic.BE.Spp_BB_TI = 10*log10(Spp_k/(data.pref^2));
%    if method.dBA == 1
%        RA      = (12200^2 * acoustic.fbb.^4) ./ ( (acoustic.fbb.^2 + 20.6^2).*sqrt((acoustic.fbb.^2 + 107.7^2).*(acoustic.fbb.^2 + 737.9^2)).*(acoustic.fbb.^2 + 12200^2) );
%        A       = 2.0 + 20*log10(RA);
%        Spp_k_A = Spp_k;
%        for  k = 1 : geom.Nb
%            Spp_k_A(k,:)        = (data.pref^2)*10.^(log10(Spp_k(k,:)/(data.pref^2))+A/10);
%        end
%        acoustic.BE.SppA_BB_TI   = 10*log10(Spp_k_A/(data.pref^2));
%    end
%end
Spp_k_A=zeros(geom.Nb,length(acoustic.fupper));
%A-weighting
if (method.freqband==0)
    for i = 1 : geom.Nb
    k=1;
    for j=1:length(acoustic.fbb)
       acoustic.BE.Spp_BB_TI(i,j) = 10*log10(Spp_k(i,j)/(data.pref^2));
       if (round(acoustic.fbb(j)) > acoustic.flower(k) && round(acoustic.fbb(j)) < acoustic.fupper(k))
         Spp_k_A(i,k)  = Spp_k_A(i,k) + (data.pref^2)*10.^( log10(Spp_k(i,j)/(data.pref^2)) );
       elseif (round(acoustic.fbb(j))==acoustic.fupper(k))
         Spp_k_A(i,k)  = 10.0*log10(Spp_k_A(i,k)/data.pref^2) + acoustic.ponderation(k);
          k=k+1
       end
    end
    OASPLA_BE_BB_TI = 0;
    for j = 1 : length(Spp_k_A(i,:)) %sum for all frequency bands (multiplied by bandwith)
        OASPLA_BE_BB_TI = OASPLA_BE_BB_TI + Spp_k_A(i,j);
    end
    acoustic.BE.OASPLA_BB_TI(i)   = 10*log10(OASPLA_BE_BB_TI/(data.pref^2));
    end
elseif (method.freqband==1)
    for i = 1 : geom.Nb
    OASPLA_BE_BB_TI = 0;
    for j = 1 : length(acoustic.fbb) %sum for all frequency bands (multiplied by bandwith)
        Spp_k_A(i,j)   = (data.pref^2)*10.^( log10(Spp_k(i,j)/(data.pref^2))+acoustic.ponderation(j)/10 );
        OASPLA_BE_BB_TI = OASPLA_BE_BB_TI + Spp_k_A(i,j);
        acoustic.BE.Spp_BB_TI(i,j) = 10*log10(Spp_k(i,j)/(data.pref^2));
        acoustic.BE.SppA_BB_TI(i,j) = 10*log10(Spp_k_A(i,j)/(data.pref^2));
    end
    acoustic.BE.OASPLA_BB_TI(i)   = 10*log10(OASPLA_BE_BB_TI/(data.pref^2));
    end
end
