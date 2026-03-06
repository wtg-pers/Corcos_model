%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Weight] = Aerodynamic_MonoAirfoil_WeightWatcher(airfoil,R,C,NBlades)
%
% CONVENTION DE SIGNE POUR LA REPRESENTATION DU PROFIL : 
% Le profil est centre et oriente de maniere a avoir 
% le bord de fuite à gauche et le bord d'attaque à droite. 
%
fid = fopen(['./AirfoilDataBase/' airfoil '.dat'],'r');
%fid = fopen(['./airfoil/',airfoil,'.dat'],'r');
M = textscan(fid, '%f %f' ,'delimiter',' \t');
fclose(fid);
Xp = M{1};
Yp = M{2};
i=1;
while (Xp(i) > min(Xp(:)))
X1(i) = Xp(i);
Y1(i) = Yp(i);
i=i+1;
end
X1(i) = Xp(i);
Y1(i) = Yp(i);
iii=0;
for ii=i+1:length(Xp)
iii=iii+1;
X2(iii) = Xp(ii);
Y2(iii) = Yp(ii);
end
NF=length(Xp);
N1=length(X1);N2=length(X2);

NP=length(C);
dR = zeros(NP,1);
dR(1) = (R(2) - R(1)) * 0.5; dR(end) = (R(end) - R(end-1)) * 0.5;
for k=2:NP-1
   dR(k) = (R(k+1) - R(k-1)) * 0.5;
end
%
for j=1:NP
   % EXTRADOS
   for i=1:N1
      X1_bis(i,j) = C(j)*X1(i);
      Y1_bis(i,j) = C(j)*Y1(i);
   end
   % INTRADOS
   for i=1:N2
      X2_bis(i,j) = C(j)*X2(i);
      Y2_bis(i,j) = C(j)*Y2(i);
   end
end

% INITIALISATION ET DECALAGE VERTICAL DU PROFIL;
Y1_bis=Y1_bis+1; Y2_bis=Y2_bis+1; S = zeros(NP,1); V = 0;
%
for k=1:NP
   S1 = 0; S2 = S1 ; 
   %
   S1 = S1 + Y1_bis(1,k)  * abs(X1_bis(2,k)  - X1_bis(1,k)   )*0.5 ; 
   S1 = S1 + Y1_bis(N1,k) * abs(X1_bis(N1,k) - X1_bis(N1-1,k))*0.5 ; 
   for j=2:N1-1
      S1 = S1 + Y1_bis(j,k) * abs(X1_bis(j+1,k) - X1_bis(j-1,k))*0.5 ; 
   end
   %
   S2 = S2 + Y2_bis(1,k)  * abs(X2_bis(2,k)  - X2_bis(1,k)   )*0.5 ; 
   S2 = S2 + Y2_bis(N2,k) * abs(X2_bis(N2,k) - X2_bis(N2-1,k))*0.5 ; 
   for j=2:N2-1
      S2 = S2 + Y2_bis(j,k) * abs(X2_bis(j+1,k) - X2_bis(j-1,k))*0.5 ; 
   end
   %
   S(k) = abs(S1-S2);
   V = V + S(k)*dR(k);
end
W = V*1000*1.10;
%
Weight = W * NBlades;
