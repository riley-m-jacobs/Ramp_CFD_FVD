function [HLLE,Se,error] = HLLE(primL, primR, n, gamma)

%
% This is the function for the HLLE
% --inputs--
% primL(1:4) = left state (rhoL, uL, vL, pL)
% primR(1:4) = right state (rhoR, uR, vR, pR)
% n          = unit normals
% gamma      = Ratio of specific heats for gas (using air)
%         ...
%
% --outputs--
% HLLE(1:4)   = numerical flux for HLLE scheme
% Se         = wave speed vector
% error      = indicator if P < 0 --> will show an error quickly!
%

nx = n(1);
ny = n(2);

error = 0;

[FLx,FLy,VL,aL,errorL] = eulerflux(primL, gamma);
[FRx,FRy,VR,aR,errorR] = eulerflux(primR, gamma);

if errorL == 1 || errorR == 1 % indicator if P < 0 --> will show an error quickly!
    error = 1;
end

uL = VL(1);
vL = VL(2);

uR = VR(1);
vR = VR(2);

unL = uL * nx + vL * ny; 
unR = uR * nx + vR * ny;

FL = FLx*nx + FLy*ny;
FR = FRx*nx + FRy*ny;

sLmin = min([0 unL - aL]); 
sRmin = min([0 unR - aR]);
sLmax = max([0 unL + aL]); 
sRmax = max([0 unR + aR]);

smin = min([sLmin sRmin]); 
smax = max([sLmax sRmax]); 

HLLE = 0.5*(FR+FL) - 0.5*((smax+smin)/(smax-smin))*(FR-FL) + ((smax*smin)/(smax-smin))*(primR-primL);

Se = max([abs(unL)+aL abs(unR)+aR]);

end

