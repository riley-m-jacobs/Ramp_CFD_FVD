function [Roe, Se, error] = Roe(primL,primR,n,gamma)

% --inputs--
% primL(1:4) = left state (rhoL, uL, vL, pL)
% primR(1:4) = right state (rhoR, uR, vR, pR)
% n          = unit normals
% gamma      = Ratio of specific heats for gas (using air)
%         ...
%
% --outputs--
% Roe(1:4)   = numerical flux for Roe scheme
% Se         = wave speed vector
% error      = indicator if P < 0 --> will show an error quickly!

primL = primL';
primR = primR';

error = 0;

% normal vectors
nx = n(1);
ny = n(2);

% Tangent vectors
mx = -ny;
my = nx;

% Left state
rL = primL(1);
uL = primL(2)/rL;
vL = primL(3)/rL;
vnL = uL*nx+vL*ny;
vtL = uL*mx+vL*my;
pL = (gamma-1)*( primL(4) - rL*(uL^2+vL^2)/2 );
%aL = sqrt(gamma*pL/rL);
HL = ( primL(4) + pL ) / rL;

% Right state
rR = primR(1);
uR = primR(2)/rR;
vR = primR(3)/rR;
vnR = uR*nx+vR*ny;
vtR = uR*mx+vR*my;
pR = (gamma-1)*( primR(4) - rR*(uR^2+vR^2)/2 );
%aR = sqrt(gamma*pR/rR);
HR = ( primR(4) + pR ) / rR;

% First compute the Roe Averages
RT = sqrt(rR/rL);
r = RT*rL;
u = (uL+RT*uR)/(1+RT);
v = (vL+RT*vR)/(1+RT);
H = ( HL+RT* HR)/(1+RT);
a = sqrt( (gamma-1)*(H-(u^2+v^2)/2) );
vn = u*nx+v*ny;
vt = u*mx+v*my;

% Wave Strengths
dr = rR - rL;     
dp = pR - pL;     
dvn= vnR - vnL;     
dvt= vtR - vtL;

dV = [(dp-r*a*dvn )/(2*a^2); 
    r*dvt/a; dr-dp/(a^2); 
    (dp+r*a*dvn)/(2*a^2)];

% Wave Speed
ws = [abs(vn-a); 
    abs(vn); 
    abs(vn); 
    abs(vn+a)];

% Harten's Entropy Fix JCP(1983), 49, pp357-393:
% only for the nonlinear fields.
dws(1)=1/5; 

if ws(1) < dws(1)
    ws(1)=(ws(1)*ws(1)/dws(1)+dws(1))/2; 
end

dws(4) = 1/5; 

if ws(4) < dws(4) 
    ws(4) = (ws(4)*ws(4)/dws(4)+dws(4))/2; 
end

% Right Eigenvectors       
Rv = [1, 0, 1,  1;
  u-a*nx, a*mx, u, u+a*nx;
  u-a*ny, a*my, u, u+a*ny;
  H-vn*a, vt*a,(u^2+v^2)/2, H+vn*a];

% Left and Right fluxes
FL=[rL*vnL; 
    rL*vnL*uL + pL*nx; 
    rL*vnL*vL + pL*ny; 
    rL*vnL*HL];

FR=[rR*vnR; 
    rR*vnR*uR + pR*nx; 
    rR*vnR*vR + pR*ny; 
    rR*vnR*HR];

% Dissipation Term
Roe = (FL + FR - Rv*(ws.*dV))/2;
Roe = Roe';

Se = max(abs(ws));

end
