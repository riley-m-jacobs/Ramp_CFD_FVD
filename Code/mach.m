function M = mach(u, gamma)

% --inputs--
% u(1:4)     = state (rhoL, uL, vL, pL)
% gamma      = Ratio of specific heats for gas (using air)
%         ...
%
% --outputs--
% M          = mach number M for the euler state u

U = 1./u(:,1).*sqrt(u(:,2).^2+u(:,3).^2); 
p = (gamma-1).*(u(:,4)-0.5.*u(:,1).*U.^2); 
a = (gamma*p./u(:,1)).^(0.5);
M = U./a; 

end

