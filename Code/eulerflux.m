function [Fx,Fy,uvec,a,error] = eulerflux(u,gamma)

% --inputs--
% u(1:4)     = state (rhoL, uL, vL, pL)
% gamma      = Ratio of specific heats for gas (using air)
%         ...
%
% --outputs--
% Fx         = flux in the x-direction
% Fy         = flux in the y-direction
% uvec       = components of velocity
% a          = speed of sound
% error      = indicator if P < 0 --> will show an error quickly!

error = 0; 

uvec = [u(2)/u(1) u(3)/u(1)];
U = (norm(uvec))^2;
P = (gamma-1)*(u(4) - 0.5*u(1)*U);

if P < 0 %throw error if pressure is negative 
    error = 1;
end

H = u(4)/u(1) + P/u(1);
a = sqrt(gamma*P/u(1));

Fx = [u(2) u(2)^2/u(1) + P (u(2)*u(3))/u(1) u(2)*H];
Fy = [u(3) (u(2)*u(3))/u(1) u(3)^2/u(1) + P u(3)*H];

end

