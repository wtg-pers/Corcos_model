%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trailing-edge noise computation 
% Model from Rozenberg et al. AIAA 2010 
% Wall pressure spectrum by ...
%
% V1.0 by R.SERRE and N.PICHON 23/09/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [acoustic] = Acoustic_broadband_TE(test,geom,flow,acoustic,obs,data,method)

% ==== Corcos spanwise coherence options (defaults) ======================
%%% >>> CORCOS UPDATE
if ~isfield(method,'CorcosFlag'),   method.CorcosFlag   = 0;   end  % 0: off, 1: on
if ~isfield(method,'CorcosAlpha'),  method.CorcosAlpha  = 0.15; end % Corcos 감쇠 계수 (~0.1~0.3)
%%% <<< CORCOS UPDATE
%
% Load Boundary layer informations
if (method.airfopt==0)
    boundary    = Aerodynamic_MonoAirfoil_BoundaryLayer(data,geom);
end
%
%Initialisation
Nb_bl  = geom.Nb;
Nfreq  = length(acoustic.fbb);

Spp         = zeros(Nfreq,1);                  % 최종(기존) TE 스펙트럼
Spp_k       = zeros(length(geom.Rb),Nfreq);    % 각 블레이드 요소별 스펙트럼

%%% >>> CORCOS UPDATE : Corcos 합산용 컨테이너
Spp_cor     = zeros(Nfreq,1);                  % Corcos 모델 TE 스펙트럼
%%% <<< CORCOS UPDATE

Doppler     = zeros(acoustic.Ntobs,1);
Spp_psi     = zeros(acoustic.Ntobs,1);
PSI         = [0:acoustic.dtobs:acoustic.Tobs]*test.RPM*2*pi/60;    %angular step
%
sinTHETA    = obs.x0(1)/obs.Robs;           % theta is the azimuthal angle
%
%calculations
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
    Radius   = geom.Rb(k);               % Rayon de l'element de pale
    Chord    = geom.Cb(k);               % Chord de l'element de pale
    U        = test.RPM*2*pi/60*Radius;  % blade element speed
    Rc       = U*Chord / data.nu;        % local Reynolds
    Mach     = U/data.c;                 % Mach number
    Mt       = Mach;                     % Rotational Mach number
    beta     = 1 - (Mach)^2;             % Prandtl-Glauert Compressibility coefficient
    L        = geom.dR;                  % Blade element span-wise length

    Uc       = 0.75*U;                   % convective speed for Rozenberg
    bc       = 1.4;
    %
    % load corresponding boundary layer values
    [~,indexOfClosestReynolds] = min(abs(boundary.Re-Rc)); % find closest available index in boundary layer tables
    [~,indexOfLocalAlpha]      = min(abs(boundary.Alpha - flow.AOA(k)));
    %
    dsp = boundary.dsp(indexOfLocalAlpha,indexOfClosestReynolds) * Chord ;
    dss = boundary.dss(indexOfLocalAlpha,indexOfClosestReynolds) * Chord ;
    Uep = boundary.Uep(indexOfLocalAlpha,indexOfClosestReynolds) / Chord ;
    Ues = boundary.Ues(indexOfLocalAlpha,indexOfClosestReynolds) / Chord ;
    Cfp = boundary.Cfp(indexOfLocalAlpha,indexOfClosestReynolds);
    Cfs = boundary.Cfs(indexOfLocalAlpha,indexOfClosestReynolds);
    %
    if (method.freqband==0)
       startf=2;
    elseif (method.freqband==1)
       startf=1;
    end
    %
    for i=startf:Nfreq % for each frequency band
        wobs = 2*pi*acoustic.fbb(i); %associated pulsation (from observer point of view)
        for j=1:acoustic.Ntobs %for each timestep
            Doppler(j) = 1 + Mt*sin(PSI(j))*sinTHETA - Mt*cos(PSI(j))*obs.x0(2)/obs.Robs; % Doppler factor due to blade rotation
            % observer's position in the moving frame
            Y1 = obs.x0(1) - geom.Rb(k)*cos(PSI(j));
            Y2 = obs.x0(2) - geom.Rb(k)*sin(PSI(j));
            Y3 = obs.x0(3);
            S0 = sqrt( Y1^2 + beta*(Y2^2 + Y3^2) );
            %
            wsource = Doppler(j)*wobs; % Source actual emission pulsation
            %
            % WALL PRESSURE SPECTRUM MODEL
            % Model de Kim-George (Airfoil Models) - Rozenberg et al.
            % intrados -------------------------------------------
            fac     = wsource * dsp / U ;
            if (fac < 0.06)
               PHI_p = ( (1.732e-3 * fac)/(1 - 5.489*fac + 36.74*fac^2 + 0.1505*fac^5) );
            else
               PHI_p = ( (1.4216e-3 * fac)/(0.3261 + 4.1837*fac + 22.818*fac^2 + 0.0013*fac^3 + 0.0028*fac^5) );
            end
            PHI_p = PHI_p*dsp*(0.5*data.rho*U^2)^2 / U;
            % extrados -------------------------------------------
            fac     = wsource * dss / U ;
            if (fac < 0.06)
               PHI_s = ( (1.732e-3 * fac)/(1 - 5.489*fac + 36.74*fac^2 + 0.1505*fac^5) );
            else
               PHI_s = ( (1.4216e-3 * fac)/(0.3261 + 4.1837*fac + 22.818*fac^2 + 0.0013*fac^3 + 0.0028*fac^5) );
            end
            PHI_s = PHI_s*dss*(0.5*data.rho*U^2)^2 / U;
            %
            PHIpp =  PHI_p + PHI_s;
            %
            % PSD model by Rozenberg et al.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ly     = bc*Uc/wsource;
            K1     = wsource/Uc;
            mu     = K1 * Mach / beta;
            fB     = K1 + mu*(Mach+1);
            fC     = K1 + mu*(Y1/S0 - Mach);
            
            [cfren1,sfren1] = fresnelMoser(sqrt(2*(fB-fC)));Es1=cfren1-1i*sfren1;
            [cfren2,sfren2] = fresnelMoser(sqrt(2*fB));Es2=cfren2-1i*sfren2;
            %
            K      = wsource / U ;
            eps    = 1/sqrt( 1 + 1/(4*mu) );
            fD     = mu*( 1-Y1/S0 );
            fT     = sqrt( (K1+mu*(Mach+1))/(K+mu*(Mach+1)) );
            fH     = ((1+1i)*exp(-4*1i*mu)*(1-fT^2)) / (2*sqrt(pi)*(U/Uc-1)*K1*sqrt(fB));
            %
            [cfren3,sfren3] = fresnelMoser(sqrt(4*mu));Es3=cfren3-1i*sfren3;
            [cfren4,sfren4] = fresnelMoser(sqrt(2*fD));Es4=cfren4-1i*sfren4;
            [cfren5,sfren5] = fresnelMoser(sqrt(4*mu));Es5=cfren5+1i*sfren5;
            %
            fG     = (1+eps)*exp( 1i*( 2*mu + fD))*( (sin(fD-2*mu))/(fD-2*mu) ) ...
                + (1-eps)*exp( 1i*(-2*mu + fD))*( (sin(fD+2*mu))/(fD+2*mu) ) ...
                + ( ((1+eps)*(1-1i))/(2*(fD-2*mu)) )*exp( 4*1i*mu )*Es3 ...
                - ( ((1-eps)*(1+1i))/(2*(fD+2*mu)) )*exp(-4*1i*mu )*Es5 ...
                + 0.5*exp(2*1i*fD)*sqrt( 2*mu/fD )*Es4 ...
                * ( ( ((1-eps)*(1+1i))/(fD+2*mu) ) - ( ((1+eps)*(1-1i))/(fD-2*mu)) ) ;
            
            % L1 : Transfer function by Amiet
            L1     = - (exp(2*1i*fC)/(1i*fC))*(...
                (1+1i)*exp(-2*1i*fC)*sqrt( fB / (fB-fC))*Es1 - (1+1i)*Es2 + 1);
            % L2 : back-scattering correction (non-finite chord length) by Roger et Moreau
            L2     = fH*(  (exp(4*eps*1i*mu)*(1-(1+eps*1i)*(Es3))) - exp(2*1i*fD) + 1i*(fD+K+mu*(Mach-1))*fG );
            
            Lcorr  = L1 + L2;
            
            Spp_psi(j) = ((wsource*Chord*Y3)/(2*pi*data.c*S0^2))^2 * (L/2) * abs(Lcorr)^2 * PHIpp * ly;
        end
        % 각 블레이드 요소 k, 주파수 i 에 대한 PSD (Nb 개의 블레이드 포함)
        Spp_k(k,i)   = Spp_k(k,i) + geom.NBlades * trapz(Doppler.*Spp_psi)/(2*pi); %time integration

        % !!! 여기서는 더 이상 Spp(i)를 누적하지 않음 (합산은 아래에서 수행)
        % 기존 코드:
        % Spp(i) = (sqrt(Spp(i)) + sqrt(abs(Spp_k(k,i))))^2;
    end
end

%%% >>> CORCOS UPDATE : 블레이드 요소들을 합치는 단계 (기존 vs Corcos)
% 기존 모델 : 모든 블레이드 요소가 완전 상관 (γ = 1)
for i = 1:Nfreq
    A = sqrt( max(Spp_k(:,i),0) );  % 각 요소의 복소 진폭 크기 (위상은 무시)
    Spp(i) = (sum(A))^2;
end

% Corcos 모델 : spanwise 거리에 따라 coherence가 감소
if method.CorcosFlag
    Nb    = geom.Nb;
    Rbvec = geom.Rb(:);  % 반경 위치 [m]
    for i = 1:Nfreq
        A = sqrt( max(Spp_k(:,i),0) );
        omega = 2*pi*acoustic.fbb(i);
        % 자기 항 (k == m)
        Spp_cor(i) = sum(A.^2);
        % 교차 항 (k < m)
        for k = 1:Nb
            for m = k+1:Nb
                dZ   = abs(Rbvec(k) - Rbvec(m));           % spanwise(반경) 거리
                Uc_k = 0.75*(test.RPM*2*pi/60*Rbvec(k));
                Uc_m = 0.75*(test.RPM*2*pi/60*Rbvec(m));
                Uc_km = 0.5*(Uc_k + Uc_m);
                gamma_km = exp( - method.CorcosAlpha * omega * dZ / Uc_km );
                Spp_cor(i) = Spp_cor(i) + 2*A(k)*A(m)*gamma_km;
            end
        end
    end
end
%%% <<< CORCOS UPDATE

%
%Overall trailing edge noise calculation
OASPL_BB_TE = 0;
OASPL_BB_TE_Corcos = 0;   %%% >>> CORCOS UPDATE

if (method.freqband==0)
    bw = ones (1,Nfreq); % bandwith is constant equals to 1Hz if tonal noise frequency spacing was selected
else
    bw = zeros(1,Nfreq);
    for i = 1 : Nfreq-1 % determines bandwidth when freq bands are equally spaced in the log space
        bw(i) = 2*(10.^((log10(acoustic.fbb(i+1))+log10(acoustic.fbb(i)))/2)-acoustic.fbb(i));
    end
    bw(Nfreq) = 2*(-10.^((log10(acoustic.fbb(end))+log10(acoustic.fbb(end-1)))/2)+acoustic.fbb(end));
end
%
Spp   = Spp.*bw';
for i = 1 : geom.Nb
    Spp_k(i,:) = Spp_k(i,:).*bw(:)';
end

if method.CorcosFlag
    Spp_cor = Spp_cor .* bw';
end

%
for i = 1 : Nfreq %sum for all frequency bands (multiplied by bandwith)
    OASPL_BB_TE = OASPL_BB_TE + Spp(i);
    if method.CorcosFlag
        OASPL_BB_TE_Corcos = OASPL_BB_TE_Corcos + Spp_cor(i);
    end
end
%OASPL_BB_TE
acoustic.OASPL_BB_TE = 10*log10(OASPL_BB_TE/(data.pref^2));
if method.CorcosFlag
    acoustic.OASPL_BB_TE_Corcos = 10*log10(OASPL_BB_TE_Corcos/(data.pref^2));
end
%
Spp_A=zeros(length(acoustic.fupper),1);
%A-weighting (기존 TE용)
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
    OASPLA_BB_TE = 0;
    for i = 1 : length(Spp_A)
        OASPLA_BB_TE = OASPLA_BB_TE + Spp_A(i);
    end
    acoustic.OASPLA_BB_TE   = 10*log10(OASPLA_BB_TE/(data.pref^2));
elseif (method.freqband==1)
    OASPLA_BB_TE = 0;
    for i = 1 : length(acoustic.fbb)
        Spp_A(i)   = (data.pref^2)*10.^( log10(Spp(i)/(data.pref^2))+acoustic.ponderation(i)/10 );
        OASPLA_BB_TE = OASPLA_BB_TE + Spp_A(i);
    end
    acoustic.OASPLA_BB_TE   = 10*log10(OASPLA_BB_TE/(data.pref^2));
end
%export data
acoustic.Spp_BB_TE  = 10*log10(Spp/(data.pref^2));
acoustic.SppA_BB_TE = 10*log10(Spp_A/(data.pref^2));

%%% >>> CORCOS UPDATE : Corcos A-weighted 결과
if method.CorcosFlag
    Spp_A_cor = zeros(size(Spp_cor));
    if (method.freqband==0)
        k = 1;
        for j=1:length(acoustic.fbb)
            if (round(acoustic.fbb(j)) > acoustic.flower(k) && round(acoustic.fbb(j)) < acoustic.fupper(k))
                Spp_A_cor(k)  = Spp_A_cor(k) + (data.pref^2)*10.^( log10(Spp_cor(j)/(data.pref^2)) );
            elseif (round(acoustic.fbb(j))==acoustic.fupper(k))
                Spp_A_cor(k)  = 10.0*log10(Spp_A_cor(k)/data.pref^2) + acoustic.ponderation(k);
                k = k+1;
            end
        end
        OASPLA_BB_TE_Corcos = 0;
        for i = 1:length(Spp_A_cor)
            OASPLA_BB_TE_Corcos = OASPLA_BB_TE_Corcos + Spp_A_cor(i);
        end
        acoustic.OASPLA_BB_TE_Corcos = 10*log10(OASPLA_BB_TE_Corcos/(data.pref^2));
        acoustic.SppA_BB_TE_Corcos   = 10*log10(Spp_A_cor/(data.pref^2));
    elseif (method.freqband==1)
        OASPLA_BB_TE_Corcos = 0;
        for i = 1:length(acoustic.fbb)
            Spp_A_cor(i)   = (data.pref^2)*10.^( log10(Spp_cor(i)/(data.pref^2))+acoustic.ponderation(i)/10 );
            OASPLA_BB_TE_Corcos = OASPLA_BB_TE_Corcos + Spp_A_cor(i);
        end
        acoustic.OASPLA_BB_TE_Corcos = 10*log10(OASPLA_BB_TE_Corcos/(data.pref^2));
        acoustic.SppA_BB_TE_Corcos   = 10*log10(Spp_A_cor/(data.pref^2));
    end
end
%%% <<< CORCOS UPDATE

%
%Blade element post-processing (기존: 원래 TE 모델 기준)
Spp_k_A=zeros(geom.Nb,length(acoustic.fupper));
%A-weighting
if (method.freqband==0)
    for i = 1 : geom.Nb
        k=1;
        for j=1:length(acoustic.fbb)
           acoustic.BE.Spp_BB_TE(i,j) = 10*log10(Spp_k(i,j)/(data.pref^2));
           if (round(acoustic.fbb(j)) > acoustic.flower(k) && round(acoustic.fbb(j)) < acoustic.fupper(k))
             Spp_k_A(i,k)  = Spp_k_A(i,k) + (data.pref^2)*10.^( log10(Spp_k(i,j)/(data.pref^2)) );
           elseif (round(acoustic.fbb(j))==acoustic.fupper(k))
             Spp_k_A(i,k)  = 10.0*log10(Spp_k_A(i,k)/data.pref^2) + acoustic.ponderation(k);
              k=k+1;
           end
        end
        OASPLA_BE_BB_TE = 0;
        for j = 1 : length(Spp_k_A(i,:))
            OASPLA_BE_BB_TE = OASPLA_BE_BB_TE + Spp_k_A(i,j);
        end
        acoustic.BE.OASPLA_BB_TE(i)   = 10*log10(OASPLA_BE_BB_TE/(data.pref^2));
    end
elseif (method.freqband==1)
    for i = 1 : geom.Nb
        OASPLA_BE_BB_TE = 0;
        for j = 1 : length(acoustic.fbb)
            Spp_k_A(i,j)   = (data.pref^2)*10.^( log10(Spp_k(i,j)/(data.pref^2))+acoustic.ponderation(j)/10 );
            OASPLA_BE_BB_TE = OASPLA_BE_BB_TE + Spp_k_A(i,j);
            acoustic.BE.Spp_BB_TE(i,j) = 10*log10(Spp_k(i,j)/(data.pref^2));
            acoustic.BE.SppA_BB_TE(i,j) = 10*log10(Spp_k_A(i,j)/(data.pref^2));
        end
        acoustic.BE.OASPLA_BB_TE(i)   = 10*log10(OASPLA_BE_BB_TE/(data.pref^2));
    end
end
