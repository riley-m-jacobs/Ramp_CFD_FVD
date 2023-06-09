function ufs = conservative(M,alpha,gamma)

% --inputs--
% M          = Free stream Mach Number M
% alpha      = An angle of attack alpha in degrees (not really used)
% gamma      = Ratio of specific heats for gas (using air)
%         ...
%
% --outputs--
% ufs(1:4)   = Calculates the state vector of [/rho, /rho*u, /rho*v, E]

ufs = zeros(1,4);

ufs(1) = 1;
ufs(2) = M*cosd(alpha);
ufs(3) = M*sind(alpha);
ufs(4) = 1/((gamma-1)*gamma) + (M^2)/(2);

end
