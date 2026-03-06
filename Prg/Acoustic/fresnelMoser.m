% From Michael Möser
% Engineering Acosutics : an introduction to noise control
% Springer, 2nd edition, p.342
%
function [cfren1,sfren1] = fresnelMoser(xarg);
%
x = abs(xarg) / sqrt(pi/2) ;
arg = pi*(x^2)/2;
s=sin(arg);
c=cos(arg);

if x>4.4

   x4=x^4;
   x3=x^3;
      x1=0.3183099 - 0.0968/x4;
      x2=0.10132   - 0.154 /x4;
      cfren1=0.5 + x1*s/x - x2*c/x3;
      sfren1=0.5 - x1*c/x - x2*s/x3;

   if xarg<0
      cfren1=-cfren1;
      sfren1=-sfren1;
   end

else  

   a0=x;
   somme=x;
   xmul=-((pi/2)^2)*(x^4);
   an=a0;
   nend=(x+1)*20;

   for n = 0:1:nend
      xnenn=(2*n+1)*(2*n+2)*(4*n+5);
      an1=an*(4*n+1)*xmul/xnenn;
      somme=somme+an1;
      an=an1;
   end

   cfren1=somme;
   a0=(pi/6)*(x^3);
   somme=a0;
   an=a0;
   nend=(x+1)*20;
   
   for n = 0:1:nend
      xnenn=(2*n+2)*(2*n+3)*(4*n+7);
      an1=an*(4*n+3)*xmul/xnenn;
      somme=somme+an1;
      an=an1;
   end

   sfren1=somme;

   if xarg<0
      cfren1=-cfren1;
      sfren1=-sfren1;
   end

end
