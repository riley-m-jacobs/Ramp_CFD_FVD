function [R, SumR, error, dt_A] = euler_solver(method,V,E,inedges,bedges,ufs,u,R,S,CFL,gamma)

% --inputs--
% method  = Riemann Flux method of choice: 'Roe', 'HLLE' or 'HLLE+'
% V       = vertices of the mesh, V(:,1) = x and V(:,2) = y
% E       = nodes of each of the triangles in the mesh
% inedges = interior edge data for the mesh as a matrix of the form [nA nB nx ny dl Eb]
% bedges  = boundary edge data for the mesh as a matrix of the form [nA nB nx ny dl Eb] 
% ufs     = free stream state
% u       = current conservative state
% R       = residual vector
% S       = wave speed vector
% e_ind   = error indicator for interior edges
% CFL     = convergence condition by Courant–Friedrichs–Lewy
% gamma   = ratio of specific heats for gas (using air)
%         ...
%
% --outputs--
% R       = residual vector
% SumR    = cumulative sum of the elements in R_list squared
% error   = indicator if P < 0 --> will show an error quickly!
% dt_A    = change in time per each single step
% Pe      = pressure along the bottom surface (MIGHT NEED TO CHANGE THIS!! -> check with literature, am i right??)

R = R.*0;             % Initialize residuals 
S = S.*0;             % Initialize wave speeds
error = 0;            % Initialize error 

% Interior edges
for e = 1:size(inedges,1)   

    if error == 1 
        break
    end
        
    % [nL nR nx ny dl EL ER]
    EL = inedges(e,6);
    ER = inedges(e,7);
    n = [inedges(e,3) inedges(e,4)];
    dl = inedges(e,5);
    uL = u(EL,:);
    uR = u(ER,:);  

    switch method
        case 'Roe',  [F_hat,Se,error] = Roe(uL,uR,n,gamma);    
        case 'HLLE', [F_hat,Se,error] = HLLE(uL,uR,n,gamma); 
        otherwise, error('Splitting method not set.');
    end

    % Residual at boundary
    R(EL,:) = R(EL,:) + F_hat*dl;
    R(ER,:) = R(ER,:) - F_hat*dl;

    % Add wave speed 
    S(EL) = S(EL) + Se*dl;
    S(ER) = S(ER) + Se*dl;

end

% Supersonic Inlet
edges = bedges{4};
for e = 1:size(edges,1)
    if error == 1 
        break
    end

    % [nL nR nx ny dl Eb]
    Eb = edges(e,6); % Need only LHS because at inlet
    n = [edges(e,3) edges(e,4)];
    dl = edges(e,5);
    ub = u(Eb,:); 
    
    switch method
        case 'Roe',  [F_hat,Se,error] = Roe(ub,ufs,n,gamma);    
        case 'HLLE', [F_hat,Se,error] = HLLE(ub,ufs,n,gamma); 
        otherwise, error('Splitting method not set.');
    end

    % Residual at boundary
    R(Eb,:) = R(Eb,:) + F_hat*dl;
    
    % Add wave speed 
    S(Eb) = S(Eb) + Se*dl;

end

% Solid inviscid walls top/bottom
for b = [1 3]
    if error == 1 %if there was an error
        break
    end
    
    edges = bedges{b};
    for e = 1:size(edges,1)
        if error == 1 %if there was an error
            break
        end
        
        % [nL nR nx ny dl Eb]
        % nA = edges(e,1);
        % nB = edges(e,2);
        Eb = edges(e,6); 
        n = [edges(e,3) edges(e,4)];
        dl = edges(e,5);
        ub = u(Eb,:); 

        %compute numerical flux
        % U_r = U - 2*(U dot nhat))

        Vel = [ub(2)/ub(1) ub(3)/ub(1)];
        Veltan = Vel - 2*dot(Vel,n).*n;

        magV2 = (norm(Veltan))^2;
        %mag_v2 = (norm(Vel))^2;

        P = (gamma-1)*(ub(4) - 0.5*ub(1)*magV2); 
        % Is this right? Yes, want to force the flow to be tangential so calc Pressures with the tangential velocity

        if P < 0 %throw error if pressure is negative 
            error = 1;
            F_hat = 0;
            Se = 0;
        else
            F_hat = [0 P*edges(e,3) P*edges(e,4) 0];
            c = sqrt(gamma*P/ub(1));
            Se = c; 
        end

         % Residual at boundary
        R(Eb,:) = R(Eb,:) + F_hat*dl;

        % Add wave speed 
        S(Eb) = S(Eb) + Se*dl;

    end
end

% Supersonic outlet
edges = bedges{2};
for e = 1:size(edges,1)
    if error == 1 %if there was an error
        break
    end
    
    % [nL nR nx ny dl Eb]
    Eb = edges(e,6);
    %Eb = edges(e,7); 
    n = [edges(e,3) edges(e,4)];
    dl = edges(e,5);
    ub = u(Eb,:); 

    % Compute numerical flux = only euler on LHS because going out
    [Fbx,Fby,Vel,c,error] = eulerflux(ub,gamma);
    F_hat = Fbx*n(1) + Fby*n(2);
    vel = dot(Vel,n);
    Se = abs(vel)+c; 

    % Residual at boundary
    R(Eb,:) = R(Eb,:) + F_hat*dl;

    % Add wave speed 
    S(Eb) = S(Eb) + Se*dl;
end

% Convergence
SumR = 0;
nR = size(R,1)*4; 
R_list = reshape(R,[nR,1]);

% Cumulative sum of the elements in R_list squared --> in num. textbook
% pseudo code (?) maybe this is the problem?
for i = 1:nR
    SumR = SumR + R_list(i)^2;
end

for i = 1:size(E,1)
    dt_A(i) = 2*CFL/S(i); % Calculate courant number for each time step, choose min !
end

dt_A = min(dt_A);
end