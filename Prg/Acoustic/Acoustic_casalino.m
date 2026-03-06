%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Casalino
%
% V4.0 by R.SERRE, C.NANA, and N.PICHON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [acoustic] = Acoustic_casalino(acoustic,test,geom,flow,obs,data,method)

geom.Tb = geom.Tb * pi/180;

%initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=acoustic.dtobs;
ddt=2*dt;
qdt=4*dt;
sdt=6*dt;
% Time delay for derivations
Nt_bis          = length(acoustic.tobs);
acoustic.pL     = zeros(1,Nt_bis); % 10 extra samples are used because of derivation initialisation
acoustic.pT     = zeros(1,Nt_bis);
delta           = atan(geom.Fxb./geom.Rb); % init sweep angle vector
acoustic.tret   = zeros(geom.Nq,Nt_bis); % initilisation of retarded time vector
%blade element init
if method.BE == 1
    BE_pL   = zeros(geom.Nb,Nt_bis);
    BE_pT   = zeros(geom.Nb,Nt_bis);
end
%
%acoustic computations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : geom.Nb % for each blade element
    %
    %init derivation variables for the blade element%%%%%%%%%%%%%%%%%%%%%%%
    Yt=zeros(3,1);
    Ytnm1=Yt;
    Ytnm2=Yt;
    Ytnm3=Yt;
    Ytnm4=Yt;
    Lt=zeros(3,1);
    Ltnm1=Lt;
    Ltnm2=Lt;
    Ltnm3=Lt;
    Ltnm4=Lt;
    vt=zeros(3,1);
    vtnm1=vt;
    vtnm2=vt;
    vtnm3=vt;
    vtnm4=vt;
    nkt=zeros(3,1);
    nktnm1=nkt;
    nktnm2=nkt;
    nktnm3=nkt;
    nktnm4=nkt;
    rt=0;
    rtnm1=rt;
    rtnm2=rt;
    rtnm3=rt;
    rtnm4=rt;
    %
    for i = 3 : Nt_bis % for each timestep over a blade revolution (+10 samples)
        gk      = atan(-flow.InducedVelocity(k)/((test.RPM*2*pi/60)*geom.Rb(k))); % calculation of inflow angle
        x       = [0 0 test.Speed*acoustic.tobs(i)]'+obs.x0;                        % Observer's position
        eta     = [0 ; geom.Rb(k) ; 0];                                             % Panel coordinates in blade fixed frame
        % Resolution of the retarded time equation through Newton-Raphson
        % method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phit    = 0.0001;   % starting approximation for the root
        eps     = 1;        % convergence evaluation value
        tol     = 10^(-6);  % stopping criterion
        total   = 100;      % max number of iterations
        kphit   = 0;        % number of iterations
        format long;
        A       = data.c^2;
        B       = -2*data.c^2*acoustic.tobs(i);
        C       = 2*(eta(1)*x(1) + eta(2)*x(2));
        D       = 2*(eta(1)*x(2) - eta(2)*x(1));
        E       = data.c^2*acoustic.tobs(i)^2 - norm(x)^2 - norm(eta)^2 + 2*eta(3)*x(3);
        while ((eps > tol) && (kphit < total))
            cosp    = cos((2*pi*test.RPM/60)*phit);
            sinp    = sin((2*pi*test.RPM/60)*phit);
            fphi    = A*phit^2 + B*phit + C*cosp + D*sinp + E;
            dfphi   = 2*A*phit + B - (2*pi*test.RPM/60)*(C*sinp + D*cosp); % the derivative at x = x_k
            phitn   = phit - fphi/dfphi;    % a new approximation for the root
            eps     = abs(phitn - phit);    %difference between previous and current estimations of the root
            phit    = phitn;                %previous estimation of the root is replaced with new one
            kphit   = kphit+1;              %increase iterations count
        end
        
        acoustic.tret(k,i)     =  phit;
        %current time step data computation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        psi     = (2*pi*test.RPM/60)*acoustic.tret(k,i) + delta(k);   % Angle between blade's frame and x-frame in retarded time
        Tpsi    = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1]; % rotation matrix between blade and absolute referential
        y       = [0 0 test.Speed*acoustic.tret(k,i)]' + Tpsi*eta;    % Source location in the absolute referential
        % Time history storing ----------------------------------------
        Ytnm4=Ytnm3;
        Ytnm3=Ytnm2;
        Ytnm2=Ytnm1;
        Ytnm1=Yt;
        Yt=y;
        % -------------------------------------------------------------
        r=x-Ytnm2;                                                    % radiation vector in the x-frame
        ri      = r/norm(r);                                          % normalized radiation vector
        Ft      = flow.Lift(k)*sin(gk) + flow.Drag(k)*cos(gk);        % loading force (x component)
        Fn      = flow.Lift(k)*cos(gk)-flow.Drag(k)*sin(gk);          % loading force (z component)
        Leta    = -[Ft; 0; Fn];                                       % Force acting on the fluid around blade element in blade fixed frame
        L       = Tpsi*Leta;                                          % Force acting on the fluid in the absolute referential
        veta    = Tpsi^-1*[0;0;test.Speed] + cross([0;0;(2*pi*test.RPM/60)],eta); % Panel absolute velocity in blade fixed frame
        v       = Tpsi*veta;                                          % Panel absloute velocity in the absolute referential
        nketa   = [sin(geom.Tb(k)); 0; cos(geom.Tb(k))];              % blade surface normal vector in blade fixed frame
        nk      = Tpsi*nketa;                                         % blade surface normal vector in the absolute referential
        % /!\ ni=nk;
        % Time history storing ----------------------------------------
        Ltnm4=Ltnm3;
        Ltnm3=Ltnm2;
        Ltnm2=Ltnm1;
        Ltnm1=Lt;
        Lt=L;
        vtnm4=vtnm3;
        vtnm3=vtnm2;
        vtnm2=vtnm1;
        vtnm1=vt;
        vt=v;
        nktnm4=nktnm3;
        nktnm3=nktnm2;
        nktnm2=nktnm1;
        nktnm1=nkt;
        nkt=nk;
        rtnm4=rtnm3;
        rtnm3=rtnm2;
        rtnm2=rtnm1;
        rtnm1=rt;
        rt=r;
        % -------------------------------------------------------------
        %
        % Time derivation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Ordre 4 : (valeurs stockées a i-2
        if (i==3)
           % Y
           dydTau(1)=(Yt(1)-Ytnm1(1))/dt;
           dydTau(2)=(Yt(2)-Ytnm1(2))/dt;
           dydTau(3)=(Yt(3)-Ytnm1(3))/dt;
           % L
           dLdTau(1)=(Lt(1)-Ltnm1(1))/dt;
           dLdTau(2)=(Lt(2)-Ltnm1(2))/dt;
           dLdTau(3)=(Lt(3)-Ltnm1(3))/dt;
           % v
           dvdTau(1)=(vt(1)-vtnm1(1))/dt;
           dvdTau(2)=(vt(2)-vtnm1(2))/dt;
           dvdTau(3)=(vt(3)-vtnm1(3))/dt;
           % ni
           dnkdTau(1)=(nkt(1)-nktnm1(1))/dt;
           dnkdTau(2)=(nkt(2)-nktnm1(2))/dt;
           dnkdTau(3)=(nkt(3)-nktnm1(3))/dt;
           % r
           drdTau=(rt-rtnm1)/dt;
        elseif ( i==4 )
           % Y
           dydTau(1)=(Yt(1)-Ytnm2(1))/ddt;
           dydTau(2)=(Yt(2)-Ytnm2(2))/ddt;
           dydTau(3)=(Yt(3)-Ytnm2(3))/ddt;
           % L 
           dLdTau(1)=(Lt(1)-Ltnm2(1))/ddt;
           dLdTau(2)=(Lt(2)-Ltnm2(2))/ddt;
           dLdTau(3)=(Lt(3)-Ltnm2(3))/ddt;
           % v
           dvdTau(1)=(vt(1)-vtnm2(1))/ddt;
           dvdTau(2)=(vt(2)-vtnm2(2))/ddt;
           dvdTau(3)=(vt(3)-vtnm2(3))/ddt;
           % ni
           dnkdTau(1)=(nkt(1)-nktnm2(1))/ddt;
           dnkdTau(2)=(nkt(2)-nktnm2(2))/ddt;
           dnkdTau(3)=(nkt(3)-nktnm2(3))/ddt;
           % r
           drdTau=(rt-rtnm2)/ddt;
        elseif ( i>4 )
           % Y
           dydTau(1) = (4/3)*(Ytnm1(1)-Ytnm3(1))/ddt + (-1/3)*(Yt(1)-Ytnm4(1))/qdt;
           dydTau(2) = (4/3)*(Ytnm1(2)-Ytnm3(2))/ddt + (-1/3)*(Yt(2)-Ytnm4(2))/qdt;
           dydTau(3) = (4/3)*(Ytnm1(3)-Ytnm3(3))/ddt + (-1/3)*(Yt(3)-Ytnm4(3))/qdt;
           % L
           dLdTau(1) = (4/3)*(Ltnm1(1)-Ltnm3(1))/ddt + (-1/3)*(Lt(1)-Ltnm4(1))/qdt;
           dLdTau(2) = (4/3)*(Ltnm1(2)-Ltnm3(2))/ddt + (-1/3)*(Lt(2)-Ltnm4(2))/qdt;
           dLdTau(3) = (4/3)*(Ltnm1(3)-Ltnm3(3))/ddt + (-1/3)*(Lt(3)-Ltnm4(3))/qdt;
           % v
           dvdTau(1) = (4/3)*(vtnm1(1)-vtnm3(1))/ddt + (-1/3)*(vt(1)-vtnm4(1))/qdt;
           dvdTau(2) = (4/3)*(vtnm1(2)-vtnm3(2))/ddt + (-1/3)*(vt(2)-vtnm4(2))/qdt;
           dvdTau(3) = (4/3)*(vtnm1(3)-vtnm3(3))/ddt + (-1/3)*(vt(3)-vtnm4(3))/qdt;
           % ni
           dnkdTau(1) = (4/3)*(nktnm1(1)-nktnm3(1))/ddt + (-1/3)*(nkt(1)-nktnm4(1))/qdt;
           dnkdTau(2) = (4/3)*(nktnm1(2)-nktnm3(2))/ddt + (-1/3)*(nkt(2)-nktnm4(2))/qdt;
           dnkdTau(3) = (4/3)*(nktnm1(3)-nktnm3(3))/ddt + (-1/3)*(nkt(3)-nktnm4(3))/qdt;
           % r
           drdTau = (4/3)*(rtnm1-rtnm3)/ddt + (-1/3)*(rt-rtnm4)/qdt;
        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %
        M        = vtnm2/data.c;                                       % Mach speed vector in absolute frame
        Mr       = dot(M,ri);                                          % Mach speed component along radiation direction
        Lr       = dot(Ltnm2,ri);                                      % Loading component along radiation direction 
        LM       = dot(Ltnm2,M);                                       % Blade element power
        vn       = dot(vtnm2,nktnm2);                                  % panel absolute velocity in absolute referential along normal to blade
        dMdTau   = dvdTau/data.c;
        dMidTau  = dot(ri,dMdTau);
        dvdTaun  = dot(dvdTau,nktnm2);
        vdnidTau = dot(vtnm2,dnkdTau);
        dLidTau  = dot(dLdTau,ri); 
        %
        % Loading noise calcualtion
        pL1                 = dLidTau/(data.c*norm(r)*(1-Mr)^2)*geom.surface(k);
        pL2                 = (Lr-LM)/(norm(r)^2*(1-Mr)^2)*geom.surface(k);
        pL3                 = Lr*(norm(r)*dMidTau+data.c*(Mr-norm(M)^2))/(data.c*norm(r)^2*(1-Mr)^3)*geom.surface(k);
        acoustic.pL(i-2)    = acoustic.pL(i-2)+1./(4*pi)*(pL1+pL2+pL3);
        if method.BE == 1
            BE_pL(k,i-2) = BE_pL(k,i-2) + 1./(4*pi)*(pL1+pL2+pL3);
        end
        %
        % Thickness noise calculation        
        pT1                 = data.rho*(dvdTaun+vdnidTau)/(norm(r)*(1-Mr)^2)*geom.surface(k);
        pT2                 = data.rho*vn*(norm(r)*dMidTau+data.c*(Mr-norm(M)^2))...
                            /(norm(r)^2*(1-Mr)^3)*geom.surface(k);
        acoustic.pT(i-2)    = acoustic.pT(i-2)+1./(4*pi)*(pT1+pT2);
        if method.BE == 1
            BE_pT(k,i-2) = BE_pT(k,i-2) + 1./(4*pi)*(pT1+pT2);
        end
    end
end
%prepare post-processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
acoustic.tobs(1:9)          = []; %remove extra samples for future post-processing
acoustic.tobs(end-5:end)    = []; %remove extra samples for future post-processing
acoustic.pT(1:9)            = [];
acoustic.pT(end-5:end)      = [];
acoustic.pL(1:9)            = [];
acoustic.pL(end-5:end)      = [];

acoustic.p_tot  = acoustic.pL + acoustic.pT; %sum both time signal to obtain total noise

acoustic.p_tot1  = [acoustic.p_tot acoustic.p_tot(2:end)];  % p_tot1 should be twice as long as p_tot to allow up to 180° phase between blades 
acoustic.pL1     = [acoustic.pL acoustic.pL(2:end)];
acoustic.pT1     = [acoustic.pT acoustic.pT(2:end)];
acoustic.tobs1   = 0:acoustic.dtobs:2*acoustic.Tobs;           % time vector over 2 periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method.BE == 1
    for i = 1 : geom.Nb
       for j = 1 : Nt_bis;
          BE_pL_R(j) = BE_pL(i,j);
          BE_pT_R(j) = BE_pT(i,j);
       end
       BE_pT_R(1:9)          = [];
       BE_pT_R(end-5:end)    = [];
       BE_pL_R(1:9)          = [];
       BE_pL_R(end-5:end)    = [];
       %
       for j=1:length(BE_pT_R)
          acoustic.BE.pL(i,j)           = BE_pL_R(j) ;
          acoustic.BE.pT(i,j)           = BE_pT_R(j) ;
       end
    end
    %
    acoustic.BE.p_tot        = acoustic.BE.pL + acoustic.BE.pT;

    acoustic.BE.pL1     = [acoustic.BE.pL(:,:)    acoustic.BE.pL(:,2:end)   ];
    acoustic.BE.pT1     = [acoustic.BE.pT(:,:)    acoustic.BE.pT(:,2:end)   ];
    acoustic.BE.p_tot1  = [acoustic.BE.p_tot(:,:) acoustic.BE.p_tot(:,2:end)];
end
