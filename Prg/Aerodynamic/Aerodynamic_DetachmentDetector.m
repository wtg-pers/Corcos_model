%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detachment_detector                                                     %
%                                                                         %
% V1.0 by N.PICHON 13/10/2016                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [eliminate] = Aerodynamic_DetachmentDetector(geom,flow,test)

load(['./boundary_layer/',geom.airfoil,'/stall_area.mat']); % load stall area

re_bound        = 5000 : 5000 : 60000;
alpha_bound     = 0 : 0.5 : 30;
k               = 1;
detached        = 1;
while k < length(geom.Rb) && detached == 1
    [~,k1]          = min(abs(re_bound - flow.LocalReynolds(k)));
    [~,k2]          = min(abs(alpha_bound -flow.AOA(k)));
    if stall_area(k1,k2) == 0
        detached    = 0;
    else
        k           = k+1;
    end
end

if k > length(geom.Rb)*test.stall_percent % if more than test.stall_percent% of blade is stalled
    eliminate       = 1;
else
    eliminate       = 0;
end

