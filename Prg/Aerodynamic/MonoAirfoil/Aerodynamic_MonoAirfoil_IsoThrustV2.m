function [flow,test,geom] = Aerodynamic_MonoAirfoil_IsoThrustV2(test,geom,polar,data,method)


%initialisation
Yobj        = sqrt(test.thrust);
Ynm1        = 0;
Xnm1        = 0;
X           = 20000;
iter        = 0;
test.RPM    = X;
eps         = test.eps;
itmax       = test.itmax;
a = 0;
b = 0;


[flow,geom] = Aerodynamic_MonoAirfoil_BEMT(test,geom,polar,data,method);

Y=sqrt(flow.Thrust);
Xnp1 = X;
Ynp1 = Y;


while ( abs(Y - Yobj) >= eps )
   iter = iter + 1;
   %
   K = (Ynp1 - Ynm1) / (Xnp1 - Xnm1) ;
   X = ( Yobj - Ynm1 + K*Xnm1 ) / K;
   %
   test.RPM=X;
   %  
   [flow,geom] = Aerodynamic_MonoAirfoil_BEMT(test,geom,polar,data,method);
   Y = sqrt(flow.Thrust);
   %  
   if ( Y < Yobj )
      Xnm1 = X;
      Ynm1 = Y;
      a = a+1;
      b = 0;
      if b>=3
          Xnp1 = Xnp1 / 2;
          test.RPM = Xnp1;
          [flow,geom] = Aerodynamic_MonoAirfoil_BEMT(test,geom,polar,data,method);
          Ynp1 = sqrt(flow.Thrust);
          a = 0;
      end  
   elseif ( Y > Yobj )
      Xnp1 = X;
      Ynp1 = Y;
      a = 0;
      b = b+1;
      if b>=3
          Xnm1 = Xnm1 / 2;
          test.RPM = Xnm1;
          [flow,geom] = Aerodynamic_MonoAirfoil_BEMT(test,geom,polar,data,method);
          Ynm1 = sqrt(flow.Thrust);
          b = 0;
      end          
   end
   % 
   if ( iter >= itmax )
      fprintf(' WARNING : Iteration limit reached. Thrust computation stopped. \n')
      fprintf(' RPM : %i Thrust : %i \n',X,Y^2)
      test.error = 1;
      Y = Yobj;
   end
   %
end
test.it = iter;
end
