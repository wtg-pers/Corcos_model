function [Cd_loc,Cl_loc]=Aerodynamic_MultiAirfoil_PolarInterpolation(localAlpha,localReynolds,polar,j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D polar interpolation
%
% V1.0 by N.PICHON 03/10/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cl      = polar.airfopt(j).Cl;
Cd      = polar.airfopt(j).Cd;
Alpha   = polar.airfopt(j).Alpha;
Re      = polar.airfopt(j).Re;

% fin closest index%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,indexOfClosestReynolds]  = min(abs(Re-localReynolds));
[~,indexOfLocalAlpha]       = min(abs(degtorad(Alpha(:)) - localAlpha));

%read closets index values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dalpha      = abs(degtorad(Alpha(indexOfLocalAlpha)) - localAlpha)/degtorad((Alpha(end)-Alpha(1))/(length(Alpha)-1));
dRe         = abs(Re(indexOfClosestReynolds) - localReynolds)/(Re(end)-Re(1))/(length(Re)-1);
al          = Cl(indexOfLocalAlpha,indexOfClosestReynolds);
ad          = Cd(indexOfLocalAlpha,indexOfClosestReynolds);

%compute distance to the 4 closest points in the grid (respectively the
%closest bounds) and make a linear ponderation to get the actual value of
%Cl and Cd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if localReynolds< Re(1) || localReynolds>Re(end) || localAlpha<degtorad(Alpha(1)) || localAlpha>degtorad(Alpha(end)) %out of boundaries of polar surface
    if localReynolds< Re(1) || localReynolds>Re(end) % out of Reynolds boundaries
        if localAlpha<degtorad(Alpha(1)) || localAlpha>degtorad(Alpha(end)) % out of Alpha and reynolds boundaries
            bl      = al;
            bd      = ad;
            cl      = al;
            cd      = ad;
            dl      = al;
            dd      = ad;
        else % out of reynolds boundaries only
            bl      = al;
            bd      = al;
            if degtorad(Alpha(indexOfLocalAlpha))>localAlpha
                cl  = Cl(indexOfLocalAlpha-1,indexOfClosestReynolds);
                cd  = Cd(indexOfLocalAlpha-1,indexOfClosestReynolds);
                dl  = cl;
                dd  = cd;
            else
                cl  = Cl(indexOfLocalAlpha+1,indexOfClosestReynolds);
                cd  = Cd(indexOfLocalAlpha+1,indexOfClosestReynolds);
                dl  = cl;
                dd  = cd;
            end
        end
    else % out of alpha bounds only
        cl      = al;
        cd      = ad;
        if Re(indexOfClosestReynolds)>localReynolds
            bl  = Cl(indexOfLocalAlpha,indexOfClosestReynolds-1);
            bd  = Cd(indexOfLocalAlpha,indexOfClosestReynolds-1);
            dl  = bl;
            dd  = bd;
        else
            bl  = Cl(indexOfLocalAlpha,indexOfClosestReynolds+1);
            bd  = Cd(indexOfLocalAlpha,indexOfClosestReynolds+1);
            dl  = bl;
            dd  = bd;
        end
        
    end
    
else   %within boundaries of polar surface
    if degtorad(Alpha(indexOfLocalAlpha))>localAlpha && Re(indexOfClosestReynolds)>localReynolds
        bl  = Cl(indexOfLocalAlpha,indexOfClosestReynolds-1);
        bd  = Cd(indexOfLocalAlpha,indexOfClosestReynolds-1);
        cl  = Cl(indexOfLocalAlpha-1,indexOfClosestReynolds-1);
        cd  = Cd(indexOfLocalAlpha-1,indexOfClosestReynolds-1);
        dl  = Cl(indexOfLocalAlpha-1,indexOfClosestReynolds);
        dd  = Cd(indexOfLocalAlpha-1,indexOfClosestReynolds);
    elseif degtorad(Alpha(indexOfLocalAlpha))>localAlpha && Re(indexOfClosestReynolds)<localReynolds
        bl  = Cl(indexOfLocalAlpha,indexOfClosestReynolds+1);
        bd  = Cd(indexOfLocalAlpha,indexOfClosestReynolds+1);
        cl  = Cl(indexOfLocalAlpha-1,indexOfClosestReynolds+1);
        cd  = Cd(indexOfLocalAlpha-1,indexOfClosestReynolds+1);
        dl  = Cl(indexOfLocalAlpha-1,indexOfClosestReynolds);
        dd  = Cd(indexOfLocalAlpha-1,indexOfClosestReynolds);
    elseif degtorad(Alpha(indexOfLocalAlpha))<localAlpha && Re(indexOfClosestReynolds)>localReynolds
        bl  = Cl(indexOfLocalAlpha,indexOfClosestReynolds-1);
        bd  = Cd(indexOfLocalAlpha,indexOfClosestReynolds-1);
        cl  = Cl(indexOfLocalAlpha+1,indexOfClosestReynolds-1);
        cd  = Cd(indexOfLocalAlpha+1,indexOfClosestReynolds-1);
        dl  = Cl(indexOfLocalAlpha+1,indexOfClosestReynolds);
        dd  = Cd(indexOfLocalAlpha+1,indexOfClosestReynolds);
    elseif degtorad(Alpha(indexOfLocalAlpha))<localAlpha && Re(indexOfClosestReynolds)<localReynolds
        bl  = Cl(indexOfLocalAlpha,indexOfClosestReynolds+1);
        bd  = Cd(indexOfLocalAlpha,indexOfClosestReynolds+1);
        cl  = Cl(indexOfLocalAlpha+1,indexOfClosestReynolds+1);
        cd  = Cd(indexOfLocalAlpha+1,indexOfClosestReynolds+1);
        dl  = Cl(indexOfLocalAlpha+1,indexOfClosestReynolds);
        dd  = Cd(indexOfLocalAlpha+1,indexOfClosestReynolds);
    else
        bl      = al;
        bd      = ad;
        cl      = al;
        cd      = ad;
        dl      = al;
        dd      = ad;
    end
end



al      = (1-dalpha)*(1-dRe)*al;
ad      = (1-dalpha)*(1-dRe)*ad;
dl      = (dalpha)*(dRe)*dl;
dd      = (dalpha)*(dRe)*dd;
bl      = (1-dalpha)*(dRe)*bl;
bd      = (1-dalpha)*(dRe)*bd;
cl      = (dalpha)*(1-dRe)*cl;
cd      = (dalpha)*(1-dRe)*cd;

Cl_loc  = (al+bl+cl+dl);
Cd_loc  = (ad+bd+cd+dd);
