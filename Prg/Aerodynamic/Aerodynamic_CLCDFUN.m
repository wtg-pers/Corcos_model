function [CL,CD] = Aerodynamic_CLCDFUN(polar,RE,A,VSO,RMU,RHO)

CL0     = polar.CL0 ;
DCLDA   = polar.DCLDA ;
CLMIN   = polar.CLMIN ;
CLMAX   = polar.CLMAX ;
CD0     = polar.CD0 ;
CD2U    = polar.CD2U ;
CD2L    = polar.CD2L ;
CLCD0   = polar.CLCD0 ;
REREF   = polar.REREF ;
REEXP   = polar.REEXP ;


W   = RE*RMU / RHO;
MSQ = W^2 / VSO^2 ;
PG = 1 / sqrt( 1 - MSQ ) ;


%% ---- CL(alpha,Mach) function
CL = ( DCLDA * A + CL0 ) * PG ;

STALL = 0 ;
if  CL > CLMAX
    STALL = 1 ;
    ACL0 = CL0 / DCLDA ;
    CL = CLMAX * cos( A - ACL0 ) ;
elseif CL < CLMIN
    STALL = 1 ;
    ACL0 = CL0 / DCLDA ;
    CL = CLMIN * cos( A - ACL0 ) ;
end

% ASTALL = ( CLMAX - CL0 ) / DCLDA ;

%% ---- set unstalled CD
if CL > CLCD0
    CD2 = CD2U ;
else
    CD2 = CD2L ;
end

%% ---- CD-scaling factor
FAC = ( RE / REREF )^REEXP ;

%% ---- CD(CL;Re) function
CLB = CL - CLCD0 ;
CD = ( CD0 + CD2 * CLB^2 ) * FAC ;

%% ---- additional CD with normal-force stall model
if STALL == 1
    ACD0 = ( CLCD0 - CL0 ) / DCLDA ;
    DCD = 2 * sin( A - ACD0 )^2 ;
    CD = CD + DCD ;
end
