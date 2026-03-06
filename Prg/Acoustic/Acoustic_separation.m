%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEPARATION
%
% separation noise psd computation as formulated by Sinibaldi and Marino
% (attributed to Nelson and Morfley)
%
% V1.0 by  N.PICHON 10/10/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flow,acoustic] = Acoustic_separation(test,geom,flow,acoustic,obs,data,method)

%Initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Spp         = zeros(length(acoustic.f_roz),1);
Spp_k       = zeros(length(geom.Rb),length(acoustic.f_roz));
Doppler     = zeros(acoustic.Ntobs,1);
Spp_psi     = zeros(acoustic.Ntobs,1);
PSI         = [0:acoustic.dtobs:acoustic.Tobs]*test.RPM*2*pi/60;    %angular step

sinTHETA    = obs.x0(1)/obs.Robs;           % theta is the azimuthal angle


% load detached area data
load('./boundary_layer/GOE265/stall_area.mat');
%calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : length(geom.Rb)             %for each blade element
    %init values relative to blade element
    Radius   = geom.Rb(k);               % Rayon de l'element de pale
    Chord    = geom.Cb(k);               % Chord de l'element de pale
    U        = test.RPM*2*pi/60*Radius;  % blade element speed
    Mach     = U/data.c;                 % Mach number
    Mt       = Mach;                     % Rotational Mach number
    
    %check if flow is seperated or not%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Re_detached     = 5000 : 5000 : 60000; %Reynolds and Angle discretization of boundary layer stall area data 
    AoA_detached    = 0 : 0.5 : 30;
    [~,k1]          = min(abs(flow.LocalReynolds(k) - Re_detached)); %closest neighbour interpolation
    [~,k2]          = min(abs(flow.AOA(k) - AoA_detached));
    if stall_area(k1,k2) ~= 0
        flow.detached(k)    = 1;
        flow.stall_area(k)  = stall_area(k1,k2);
    else
        flow.detached(k) = 0;
    end
    %if is compute seperated noise, else set it to zero%%%%%%%%%%%%%%%%%%%%
    if flow.detached(k) == 1
        for i=1:length(acoustic.f_roz) % for each frequency band
            wobs = 2*pi*acoustic.f_roz(i); %associated pulsation (from observer point of view
            for j=1:acoustic.Ntobs %for each timestep
                if method.optim == -1
                    clc;
                    fprintf('\n%s\n','Sinibaldi and Marino separation noise calculations');
                    fprintf('\n%s\n',['blade element ',num2str(k),'/',num2str(length(geom.Rb))]);
                    fprintf('\n%s\n',['frequency band ',num2str(i),'/',num2str(length(acoustic.f_roz))]);
                    fprintf('\n%s\n',['timestep ',num2str(j),'/',num2str(acoustic.Ntobs)]);
                end
                Doppler(j) = 1 + Mt*sin(PSI(j))*sinTHETA - Mt*cos(PSI(j))*obs.x0(2)/obs.Robs; % Doppler factor due to blade rotation
                % observer's position in the moving frame
                Y1 = obs.x0(1) - geom.Rb(k)*cos(PSI(j));
                Y2 = obs.x0(2) - geom.Rb(k)*sin(PSI(j));
                Y3 = obs.x0(3);
                %
                wsource = Doppler(j)*wobs; % Source actual emission pulsation
                
                %Sinibaldi and Marino's equation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Spp_psi(j) = (wsource/(4*pi*data.c*geom.Rb(k)))^2 *(data.rho^2 * Chord * U^3*(flow.stall_area(k)*Chord*geom.dR)^2)*(Y3/geom.Rb(k))^2 ...
                    *(flow.Cd_local(k)^2/4)*8.6*10^(-7)*(60*U/(test.RPM*Chord*(1-Mt)))^3;
            end
            Spp_k(k,i)   = Spp_k(k,i) + geom.NBlades * trapz(Doppler.*Spp_psi)/(2*pi); %time integration
            
            Spp(i) = (sqrt(Spp(i)) + sqrt(abs(Spp_k(k,i))))^2; % sum over blade elements
        end
    else
        Spp_k(k,:) = 0;
    end
end

%Overall trailing edge noise calculation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OASPL_sep = 0;
if test.freq_nb == 0
    bw = ones (1,length(acoustic.f_roz)); % bandwith is constant equals to 1Hz if tonal noise frequency spacing was selected
else
    bw = zeros(1,length(acoustic.f_roz));
    for i = 1 : length(acoustic.f_roz)-1 % determines bandwidth when freq bands are equally spaced in the log space
        bw(i) = 2*(10.^((log10(acoustic.f_roz(i+1))+log10(acoustic.f_roz(i)))/2)-acoustic.f_roz(i));
    end
    bw(length(acoustic.f_roz)) = 2*(-10.^((log10(acoustic.f_roz(end))+log10(acoustic.f_roz(end-1)))/2)+acoustic.f_roz(end));
end

for i = 1 : length(acoustic.f_roz) %sum for all frequency bands (multiplied by bandwith)
    OASPL_sep = OASPL_sep + Spp(i)*bw(i);
end
acoustic.OASPL_sep = 10*log10(OASPL_sep/(data.pref^2));

%A ponderation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method.dBA == 1
    
    RA      = (12200^2 * acoustic.f_roz.^4) ./ ( (acoustic.f_roz.^2 + 20.6^2).*sqrt((acoustic.f_roz.^2 + 107.7^2).*(acoustic.f_roz.^2 + 737.9^2)).*(acoustic.f_roz.^2 + 12200^2) );
    A       = 2.0 + 20*log10(RA);
    Spp_A   = (data.pref^2)*10.^(log10(Spp/(data.pref^2))+A'/10);
    
    OASPLA_sep = 0;
    for i = 1 : length(acoustic.f_roz) %sum for all frequency bands (multiplied by bandwith)
        OASPLA_sep = OASPLA_sep + Spp_A(i)*bw(i);
    end
    acoustic.OASPLA_sep     = 10*log10(OASPL_sep/(data.pref^2));
    acoustic.Spp_A_sep      = 10*log10(Spp_A/(data.pref^2));
end
%export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acoustic.Spp_sep = 10*log10(Spp/(data.pref^2));


%Blade element post-processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method.BE == 1
    acoustic.BE.Spp_sep = 10*log10(Spp_k/(data.pref^2));
    if method.dBA == 1
        RA      = (12200^2 * acoustic.f_roz.^4) ./ ( (acoustic.f_roz.^2 + 20.6^2).*sqrt((acoustic.f_roz.^2 + 107.7^2).*(acoustic.f_roz.^2 + 737.9^2)).*(acoustic.f_roz.^2 + 12200^2) );
        A       = 2.0 + 20*log10(RA);
        Spp_k_A = Spp_k;
        for  k = 1 : geom.Nb
            Spp_k_A(k,:)        = (data.pref^2)*10.^(log10(Spp_k(k,:)/(data.pref^2))+A/10);
        end
        acoustic.BE.Spp_A_sep   = 10*log10(Spp_k_A/(data.pref^2));
    end
end
