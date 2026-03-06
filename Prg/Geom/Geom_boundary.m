function [Rb,Cb,Tb,Fxb,N] = Geom_boundary (geom) 

R   = geom.Rq;
C   = geom.C;
T   = geom.T*2*pi/360;
Fx  = geom.Fx;
dR  = geom.dR;


Rb  = R(1)+dR/2:dR:R(end)-dR/2; %elements center radial position


Cb  = spline(R,C,Rb);  %elements chord at center
Tb  = spline(R,T,Rb);  %elements twist at center
Fxb = spline(R,Fx,Rb); %elements sweep X at center
N   = length(Rb);      %number of elements center

