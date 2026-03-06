%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency spectrum 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acoustic,Lp,OASPL,LpA,OASPLA]  = Acoustic_spectrum(p,data,acoustic,test,method)
%
df               = 1/acoustic.Tsound;
f                = 0 : df : acoustic.Fs/2;      % frequency vector
acoustic.f       = f/f(end)*acoustic.Fs/2;
Lp               = zeros(length(acoustic.f),1); %init spectrum

% FFT computation
if method.window == 1    % if a different windowing is required, apply window
    p            = p.*hann(length(p))';
end

p_spectrum   = fft(p)/acoustic.Fs;  % compute spectrum

for j = 1 : length(acoustic.f)          % Sound Pressure Level
    Lp(j)    = 20*log10(abs(p_spectrum(j))/data.pref );
end
%
%OverAll Sound Pressure Level Calculation
OASPL       = 0;
for j = 2 : length(acoustic.f)
    OASPL   = OASPL  + 10^(Lp(j)*0.1)*df;
end
OASPL       = 10.0*log10(OASPL);
%
% Third Octave and A-weighting 
%
fcentre  = round( 10^3 * (2 .^ ([-18:13]'/3)) );fupper=zeros(length(fcentre),1);flower=fupper;LpA=fupper;
fd = 2^(1/6);
for j=1:length(fcentre)
   fupper(j) = round( fcentre(j) * fd );
   flower(j) = round( fcentre(j) / fd );
   RA(j) = (12200^2 * fcentre(j)^4) / ( (fcentre(j)^2 + 20.6^2)*sqrt((fcentre(j)^2 + 107.7^2)*(fcentre(j)^2 + 737.9^2))*(fcentre(j)^2 + 12200^2) );
   A(j) = 2.0 + 20*log10(RA(j));
end
%
acoustic.foct=fcentre';
acoustic.flower=flower;
acoustic.fupper=fupper;
%
if (method.dBA==0)
   acoustic.ponderation=zeros(length(A),1);
elseif (method.dBA==1)
   acoustic.ponderation=A;
end
%
%k=1;
%for j=1:length(acoustic.f)
%   if (round(acoustic.f(j)) >= flower(k) && round(acoustic.f(j)) < fupper(k))
%      LpA(k) = LpA(k) + 10^( Lp(j)*0.1 );
%   elseif (round(acoustic.f(j))>=fupper(k))
%   %elseif (round(acoustic.f(j))==fupper(k))
%      LpA(k) = 10.0*log10(LpA(k)) + acoustic.ponderation(k);
%      k=k+1;
%   end
%end
%
j=1;
k=1;
while j < length(acoustic.f)
%fprintf(' k = %i j = %i \n',k,j)
%fprintf(' Flow %5.1f F %5.1f Fhigh %5.1f \n',flower(k),acoustic.f(j),fupper(k))
   if (fupper(k)>acoustic.f(end) )
      break
   end
   if (k>=length(fupper) )
      break
   end
   if (acoustic.f(j) >= flower(k) && acoustic.f(j) < fupper(k))
      LpA(k) = LpA(k) + 10^( 0.1*Lp(j));
      j=j+1;
   elseif (acoustic.f(j)<=flower(k))
      j=j+1;
   elseif (acoustic.f(j)>=fupper(k))
      if (LpA(k)~=0)
         LpA(k) = 10.0*log10(LpA(k)) +  acoustic.ponderation(k);
         k=k+1;
      else 
         k=k+1;
      end
   end
end
%
OASPLA=0;
for j=1:length(LpA)
   OASPLA = OASPLA + 10^( LpA(j)*0.1 )*df;
end
OASPLA = 10.0*log10( OASPLA );
