function [sol, err] = driver(opts)

fig = 0;

[M2,M3] = downstream(opts.Mfs,opts.thetad,opts.gamma); 

if ( opts.sanity )
   % Analytical Sanity Check
   fprintf('M2 = %4.5f, \nM3 = %4.5f\n', M2, M3);
end

% Build the mesh
[E,V,~,n_y,bedges,inedges,x_ramp_start,x_length] = wedge_mesh(opts.n_x,opts.discretize,opts.thetad);


% Maybe start from converged solution --> Debug technique (Bisetti)
ufs = conservative(opts.Mfs,opts.alphafs,opts.gamma);

if ( opts.state )
    restartU = input('Enter restart state file name: ','s');
    u0 = fscanf(fopen(restartU, 'r'),'%f',[4,size(E,1)])';
else 
    u0 = ones(size(E,1),4).*(conservative(M2,opts.alphafs,opts.gamma));
end

% Define time step
tend = 1e6;       % Time out time for steady state
count_max = 100;  % Counter parameters for periodic state output
CFLTuner = 0;     % Initialize CFL Tuner in the off condition
    
% Discretize time domain
c0 = sqrt(opts.gamma*287*300);
vn = sqrt((opts.Mfs*c0).^2); 
lambda1 = vn + c0; 
lambda2 = vn - c0; 
CFL = 0.5;
a0 = max(abs([lambda1(:);lambda2(:)])); 
dt= CFL*min(n_y./a0, n_y./a0);

% Create initial state
Niter = 1e6;
q = u0;
Fx = 0;
Fy = 0;

% Preallocate Space
R = zeros(size(E,1),4);            % Residual vector
S = zeros(size(E,1),1);            % Wave speed vector
Rnorm = zeros(Niter,1);            % Log of residual norm 

% Solver
fprintf('\n<<<<< RAMP SOLVER START >>>>>\n');

% Initiate time trackers
count = 0;
t = 0;
N = 1;

tic
while t < tend

    if t+dt > tend
        dt=tend-t; 
    end

    t = t + dt;
    count = count + 1; % Add value to sanity check counter

    qo = q;
    
    %
    [R, SumR, error, dt_A] = euler_solver(opts.method,V,E,inedges,bedges,ufs,q,R,S,CFL,opts.gamma);
    %

    Rnorm(N) = sqrt(SumR);

    if  Rnorm(N) < 10^(-5)
        break
    end
    
    u_safe = q;

    % Check CFL for errors
    if error == 1 
        fprintf('Error occured!!! Back off on CFL! \n')
        CFLTuner = 1; % Back off on the CFL
        q = u_safe;   % Reset to previous safe state
        continue 
    end  

    CFL = opts.CFLmax - (0.5*CFLTuner);  % maybe this helps??

    q = qo - dt_A.*R; 

    % Possible 2nd Stage: (maybe help error problem??)
    %[R, SumR, error, dt_A] = euler_solver(method,V,E,inedges,bedges,ufs,q,R,S,e_ind,CFL,gamma); 

    %q = 0.75*qo + 0.25*(q - dt_A*R);

    dt = dt_A;
    
    if count >= count_max
        count = 0;
        % Output the current counter: aka N = time passed and t being the
        % current
        fprintf('Working: (N = %d & t = %d & dt = %d & dt_A = %d)\n', N, t, dt, dt_A);
        
        CFLTuner = 0; % Reset tuner
    end  

    N = N + 1;
end

if N == Niter
    fprintf('/n...DIDNT CONVERGE...\n\n');
end

fprintf('/n...time is \n\n');
cputime = toc; disp(['CPU time: ' num2str(cputime),' s'])

% Calculate converged flow properties
r = q(:,1);
U = 1./q(:,1).*sqrt(q(:,2).^2+q(:,3).^2); 
P = (opts.gamma-1).*(q(:,4)-0.5.*q(:,1).*U.^2); 
a = (opts.gamma*P./q(:,1)).^(0.5);
M = U./a; 
e = P./((opts.gamma-1)*q(:,1));   
ss = log(P./r.^opts.gamma); 

fprintf('<<<<< RAMP SOLVER COMPLETE >>>>>\n\n');

%
sol.rho = r;
sol.U_mag = U;
sol.P = P;
sol.a = a;
sol.M = M;
sol.e = e;
sol.ss = ss;
%
sol.n_x = opts.n_x;
sol.discretize = opts.discretize;
sol.n_y = n_y;
%

err = error;

if ( opts.data_save )
    outputstate = fopen('%s.txt',opts.datasaver,'w');
    fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E\n',q');
    fclose(outputstate);    
end

if ( opts.mach_plot )
    fig = fig + 1;
    figure(fig)
    clf(fig)
    meshplot(M,E,V,'Local Mach Number','Converged Mach Field',1,0);
    hold on 
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    if ( opts.save_plots )
        saveas(gcf,strcat(saveprefix,'%s Converged Mach Contours.jpg',opts.header));
    end    
end

if ( opts.pres_plot )
    fig = fig + 1;
    figure(fig)
    clf(fig)
    meshplot(P,E,V,'Local Pressure','Converged Pressure Field',1,0)
    hold on 
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    if ( opts.save_plots )
        saveas(gcf,strcat(saveprefix,'%s Converged Pressure Contours.jpg',opts.header));
    end
end

if ( opts.zoom_plot )
    fig = fig + 1;
    figure(fig)
    clf(fig)
    meshplot(M,E,V,'Local Mach Number','Converged Mach Field',1,0);
    xlim([0, (x_ramp_start+x_length)])
    ylim([0, 1])
    hold on 
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    if ( opts.save_plots )
        saveas(gcf,strcat(saveprefix,'%s Converged Mach Zoom.jpg', opts.header));
    end
end

if ( opts.mesh_plot )
    fig = fig + 1;
    figure(fig)
    clf(fig)
    hold on
    nulldat = zeros(size(E,1));
    meshplot(nulldat,E,V,'colorlabel','Full Mesh',1,1)
    colorbar('off');
    if (opts.save_plots )
        saveas(gcf,strcat(saveprefix,'Mesh Full.jpg'))
    end
end  

if ( opts.pres_bot )

    prop_in = 'Pressure';
    prop = P;
    
    % some hardcoding based on the pressure ratios of compressible flow
    % calculator
    if opts.thetad == 16.6992
        P2_P1 = 6.68634592;
    elseif opts.thetad == 30.9638
        P2_P1 = 3.12303510;
    else
        P2_P1 = 2.19465313;
    end

    [extract_line] = bottom_plot(n_y, V, E, prop);
    
    fig = fig + 1;
    figure(fig)
    plot([0,x_ramp_start], [(1/1.4),(1/1.4)],'black','linewidth',2,'DisplayName','Exact')
    hold on
    plot([x_ramp_start,x_ramp_start], [(1/1.4),(1/1.4)*P2_P1],'black','linewidth',2,'HandleVisibility','off') 
    plot([x_ramp_start+ x_length,x_ramp_start+ x_length], [(1/1.4),(1/1.4)*P2_P1],'black','linewidth',2,'HandleVisibility','off') 
    plot([x_ramp_start,x_ramp_start+ x_length], [(1/1.4)*P2_P1,(1/1.4)*P2_P1],'black','linewidth',2,'HandleVisibility','off') 
    plot([x_ramp_start+ x_length,1.2], [(1/1.4),(1/1.4)],'black','linewidth',2,'HandleVisibility','off')
    plot(extract_line(:,1),extract_line(:,3),'--','linewidth',1.5,'DisplayName','HLLE')
    legend
    ylim([0.6, 2.4])
    xlim([0.4 1.2])
    ylabel(sprintf('%s',prop_in))
    xlabel('Bottom Border')
    title(sprintf('%s along the Bottom Boundary',prop_in))

end

if ( opts.mach_bot )

    prop_in = 'Mach';
    prop = M;
    
    [extract_line] = bottom_plot(n_y, V, E, prop);

    fig = fig + 1;
    figure(fig)
    plot([0,x_ramp_start], [opts.Mfs,opts.Mfs],'black','linewidth',2,'DisplayName','Exact')
    hold on
    plot([x_ramp_start,x_ramp_start], [opts.Mfs,M2],'black','linewidth',2,'HandleVisibility','off') 
    plot([x_ramp_start,x_ramp_start+ x_length], [M2,M2],'black','linewidth',2,'HandleVisibility','off') 
    plot([x_ramp_start+ x_length,x_ramp_start+ x_length], [M3,M2],'black','linewidth',2,'HandleVisibility','off') 
    plot([x_ramp_start+ x_length,1.2], [M3,M3],'black','linewidth',2,'HandleVisibility','off')
    plot(extract_line(:,1),extract_line(:,3),'--','linewidth',1.5,'DisplayName','Roe')
    legend
    xlim([0.4 1.2])
    ylabel(sprintf('%s',prop_in))
    xlabel('Bottom Border')
    title(sprintf('%s along the Bottom Boundary',prop_in))
end

return