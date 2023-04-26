function [M2, M3] = downstream(M1,theta,gamma)
% Call the Function
[Oblique1] = rampshocks(theta,gamma,M1);

%fprintf('Results:\n')
M2 = Oblique1(4);
%fprintf("Mach number after oblique shock: %.4f\n", M2);

% Input values for P-M (NEED TO HAVE THE AEROSPACE TOOLBOX INSTALLED)
[~,nu,~] = flowprandtlmeyer(gamma,M2,'mach');
vm3 = nu + theta;
[M3,~,~] = flowprandtlmeyer(gamma,vm3,'nu');
% Display results
%fprintf("Mach number after expansion fan: %.4f\n", M3);
end

function [Oblique1] = rampshocks(theta,gamma,M1)

g = gamma;

% First Oblique 1
% Solve for beta
B_eqn = @(B) (2*cotd(B) * (M1^2*(sind(B)^2)-1)/(M1^2*(g +cosd(2*B))+2)) - tand(theta);
options = optimset('Display','off');
B0 = 20; % initial guess
beta = fsolve(B_eqn,B0,options);

Mn1      = M1*sind(beta);
Mn2      = M2(g,Mn1);
Mach2    = Mn2/sind(beta-theta);
P2P1     = 1+2*g/(g+1)*(Mn1^2-1);
P02P01   = P0(g,Mn1)/P0(g,M2(g,Mn1))*P2P1;
rho2rho1 = R0(g,M2(g,Mn1))/R0(g,Mn1)*P02P01;
T2T1     = T0(g,M2(g,Mn1))/T0(g,Mn1);

Oblique1 = [theta beta M1 Mach2 P2P1 P02P01 rho2rho1 T2T1];


function [M2_Out]  = M2(g,M1)
    M2_Out = sqrt((1+0.5*(g-1)*M1*M1)/(g*M1*M1-0.5*(g-1)));
end

function [p0] = P0(g,M)
    p0 = (1+(g-1)/2*(M^2))^(-g/(g-1));
end

function [r0] = R0(g,M)
    r0 = (1+(g-1)/2*(M^2))^(-1/(g-1));
end

function [t0] = T0(g,M)
    t0 = (1+(g-1)/2*(M^2))^(-1);
end
end

