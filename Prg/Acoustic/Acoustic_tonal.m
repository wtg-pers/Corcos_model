%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes the noise induced by the rotation
% of blades through a Ffowcs-Williams and
% Hawking analogy, using Farassat and Casalino's resolution scheme.
%
% Biblio:
% -Farassat AIAA 1981 (for FWH formulation)
% -Casalino JSV 2003, eq. (17) and (19) p. 589
% -Brentner AIAAp 1996 (for source-dominant formulation)
% -Cathy and Kelly NASA 1997 (details of computation)
% -Tam JCP 1993 (for the DRP derivation schemes)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obs,acoustic] = Acoustic_tonal(test,geom,flow,obs,data,method)
%
% Observer position
if (method.acousticTest==0)
obs.phi=obs.phi*pi/180;
obs.theta=obs.theta*pi/180;
obs.x0      =   obs.Robs*[sin(obs.phi)*cos(obs.theta) sin(obs.phi)*sin(obs.theta) cos(obs.phi)]';
elseif (method.acousticTest==1)
Xm=dlmread(['./MicrophonePosition/MICRO_' num2str(obs.micro) '_y1.dat']);
Ym=dlmread(['./MicrophonePosition/MICRO_' num2str(obs.micro) '_y2.dat']);
Zm=dlmread(['./MicrophonePosition/MICRO_' num2str(obs.micro) '_y3.dat']);
obs.x0 = [Xm Ym Zm]';
end
%
%time discretization options
timedelay = 15; 
if method.time_refinement == 0
    acoustic.Tobs    = test.cycle/(test.RPM)*60;          % Time window
    acoustic.Fsobs   = obs.Nh*(test.RPM*geom.NBlades/60); % max desired frequency to get Nh harmonics
    acoustic.Fs      = acoustic.Fsobs*2;                  % Sampling Frequency
    acoustic.Ntobs   = floor(acoustic.Tobs * acoustic.Fs);% number of samples
    acoustic.dtobs   = acoustic.Tobs/(acoustic.Ntobs-1);  % Time step
    acoustic.tobs    = 0 : acoustic.dtobs : acoustic.Tobs+timedelay*acoustic.dtobs; %10 extra samples are used to take into acount
                                                          %the discrete derivation
else
    acoustic.Tobs    = test.cycle/(test.RPM)*60;          % Time window
    acoustic.Ntobs   = obs.nppw;                          % number of samples
    acoustic.Fs      = acoustic.Ntobs/acoustic.Tobs;      % Sampling Frequency
    acoustic.Fsobs   = acoustic.Fs /2;                    % max resolved frequency
    acoustic.dtobs   = acoustic.Tobs/(acoustic.Ntobs-1);  % Time step
    acoustic.tobs    = 0 : acoustic.dtobs : acoustic.Tobs+timedelay*acoustic.dtobs; %10 extra samples are used to take into acount
                                                          %the discrete derivation
end
%
geom.surface     = geom.dR*geom.Cb;   % surface of a blade element
%
% FWH computation from Casalino's formulation
[acoustic]       = Acoustic_casalino(acoustic,test,geom,flow,obs,data,method);
%
%time signal post-processing
%
%acoustic.Ntobs
% Interpolation for the Bnumb blades
acoustic.p_tot          = zeros(1,acoustic.Ntobs);
acoustic.p_loading      = zeros(1,acoustic.Ntobs);
acoustic.p_thickness    = zeros(1,acoustic.Ntobs);
%
for i = 1 : geom.NBlades
    tint                    = acoustic.tobs + (i-1)/(test.RPM/60)/geom.NBlades; % offset time vector according to blade repartition
    acoustic.p_tot          = acoustic.p_tot + interp1(acoustic.tobs1,acoustic.p_tot1,tint,'spline'); % offset signal to first blade signal
    acoustic.p_loading      = acoustic.p_loading + interp1(acoustic.tobs1,acoustic.pL1,tint,'spline');
    acoustic.p_thickness    = acoustic.p_thickness + interp1(acoustic.tobs1,acoustic.pT1,tint,'spline');
end
%
clear acoustic.pT1; clear acoustic.pL1; clear acoustic.p_tot1; % clear temporary variables
%
%to avoid windowing problems as much as possible, if a simple fft is
%computed, tonal being periodic, many revolutions are considered (1 second
%of rotation)
acoustic.Tsound             = 1;                                        % overall sample time is 1 second
isound                      = ceil(acoustic.Tsound/acoustic.Tobs);      % number of rotor revolutions
acoustic.p_Tonal            = repmat(acoustic.p_tot(1:end-1),1,isound);   % replication of time signal over a period
acoustic.p_loading          = repmat(acoustic.p_loading(1:end-1),1,isound);
acoustic.p_thickness        = repmat(acoustic.p_thickness(1:end-1),1,isound);
%
acoustic.p_Tonal(end+1)     = acoustic.p_tot(1);
acoustic.p_loading(end+1)   = acoustic.p_loading(1);
acoustic.p_thickness(end+1) = acoustic.p_thickness(1);
%
%spectrum calculation
[acoustic,acoustic.Lp_thickness,acoustic.OASPL_thickness,acoustic.LpA_thickness,OASPLA_thickness]  = Acoustic_spectrum(acoustic.p_thickness,data,acoustic,test,method);
[acoustic,acoustic.Lp_loading,acoustic.OASPL_loading,acoustic.LpA_loading,OASPLA_loading]          = Acoustic_spectrum(acoustic.p_loading,data,acoustic,test,method);
[acoustic,acoustic.Lp_Tonal,acoustic.OASPL_Tonal,acoustic.LpA_Tonal,acoustic.OASPLA_Tonal]         = Acoustic_spectrum(acoustic.p_Tonal,data,acoustic,test,method);
%
%BE time signal post-processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method.BE == 1
    acoustic.BE.p_tot          = zeros(geom.Nb,acoustic.Ntobs);
    acoustic.BE.p_l            = zeros(geom.Nb,acoustic.Ntobs);
    acoustic.BE.p_t            = zeros(geom.Nb,acoustic.Ntobs);
    %
    for i = 1 : geom.NBlades
        for j = 1 : geom.Nb
            tint                            = acoustic.tobs + (i-1)*2*pi/geom.NBlades/(2*pi*test.RPM/60); % offset time vector according to blade repartition
%size(acoustic.tobs1)
%size(acoustic.BE.p_tot1)
%size(tint)
%pause
            acoustic.BE.p_tot(j,:)          = acoustic.BE.p_tot(j,:) + interp1(acoustic.tobs1,acoustic.BE.p_tot1(j,:),tint,'spline'); % offset signal to first blade signal
            acoustic.BE.p_l(j,:)            = acoustic.BE.p_l(j,:)   + interp1(acoustic.tobs1,acoustic.BE.pL1(j,:),tint,'spline');
            acoustic.BE.p_t(j,:)            = acoustic.BE.p_t(j,:)   + interp1(acoustic.tobs1,acoustic.BE.pT1(j,:),tint,'spline');
        end
    end
    %
    clear acoustic.BE.pT1; clear acoustic.BE.pL1; clear acoustic.BE.p_tot1; % clear temporary variables
    %
    % replication of time signal over a period
    acoustic.BE.p_Tonal              = repmat(acoustic.BE.p_tot(:,1:end),1,isound);
    acoustic.BE.p_loading            = repmat(acoustic.BE.p_l(:,1:end),1,isound);
    acoustic.BE.p_thickness          = repmat(acoustic.BE.p_t(:,1:end),1,isound);
    %
    for j = 1 : geom.Nb;
        acoustic.BE.p_Tonal(j,end+1)           = acoustic.BE.p_Tonal(j,1);
        acoustic.BE.p_loading(j,end+1)         = acoustic.BE.p_loading(j,1);
        acoustic.BE.p_thickness(j,end+1)       = acoustic.BE.p_thickness(j,1);
    end
    %
    clear acoustic.BE.p_l; clear acoustic.BE.p_t; clear acoustic.BE.p_tot; % clear temporary variables
    %spectrum calculation
    for j = 1 : geom.Nb;
        [acoustic,acoustic.BE.Lp_thickness(:,j),acoustic.BE.OASPL_thickness(j),acoustic.BE.LpA_thickness(:,j),acoustic.BE.OASPLA_thickness(j)] = Acoustic_spectrum(acoustic.BE.p_thickness(j,:),data,acoustic,test,method);
        [acoustic,acoustic.BE.Lp_loading(:,j),acoustic.BE.OASPL_loading(j),acoustic.BE.LpA_loading(:,j),acoustic.BE.OASPLA_loading(j)]         = Acoustic_spectrum(acoustic.BE.p_loading(j,:),data,acoustic,test,method);
        [acoustic,acoustic.BE.Lp_Tonal(:,j),acoustic.BE.OASPL_Tonal(j),acoustic.BE.LpA_Tonal(:,j),acoustic.BE.OASPLA_Tonal(j)]        = Acoustic_spectrum(acoustic.BE.p_Tonal(j,:),data,acoustic,test,method);
    end
end
