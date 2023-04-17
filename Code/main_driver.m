%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Euler Inviscid Compressible Flow Solver
% Compression Ramp Mesh w/ Roe, HLLE Flux Formulas
%
% Riley Jacobs
% 03.19.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close 
clc

load('colour.mat');

% Input values
Mfs = input("Enter Mach number before shock: ");
thetad = input("Enter turning angle of ramp in degrees: ");
alphafs = input("Enter direction of uniform flow: ");
gamma = input("Enter specific heat ratio: ");

% Mesh 
n_x = input("Enter the n_x or how discretize you want your mesh to be: "); % Usually unrefined is n_x = 10 and ramp is also 10
discretize = input("Define the number of nodes up the ramp: ");

% Choose method
method = input('Enter flux method: ','s'); % ROE or HLLE

% Analytical Sanity Check
[M2,M3] = downstream(Mfs,thetad,gamma); 

% Build the mesh
[E,V,B,n_y,bedges,inedges,x_ramp_start,x_length] = wedge_mesh(n_x,discretize,thetad);



%%
% Maybe start from converged solution --> Debug technique (Bisetti)
ufs = conservative(Mfs,alphafs,gamma);

answer = input('\nStart from a previously converged state? (y/n): ','s');
if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes')) 
    restartU = input('Enter restart state file name: ','s');
    u0 = fscanf(fopen(restartU, 'r'),'%f',[4,size(E,1)])';
else 
    u0 = ones(size(E,1),4).*(conservative(M2,alphafs,gamma));
end

%%
fprintf('\n');
fprintf('          IC:     \n');
fprintf('density : %2.4f \n',ufs(1));
fprintf('  x-vel : %2.4f \n',ufs(2));
fprintf('  y-vel : %2.4f \n',ufs(3));
fprintf(' Energy : %2.4f \n',ufs(4));
fprintf('\n');

% Define time step
tend = 1e6;       % Time out time for steady state
count_max = 100;  % Counter parameters for periodic state output
CFLmax = .4;      % Maximum CFL allowed for time stepping
CFLTuner = 0;     % Initialize CFL Tuner in the off condition
    
% Discretize time domain
c0 = sqrt(gamma*287*300);
vn = sqrt((Mfs*c0).^2); 
lambda1=vn+c0; 
lambda2=vn-c0; 
CFL = 0.5;
a0 = max(abs([lambda1(:);lambda2(:)])); 
dt=CFL*min(n_y./a0,n_y./a0);

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
fprintf('\n<<<<< EULER SOLVER START >>>>>\n');

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
    count = count + 1; %add value to counter

    qo = q;
    
    %
    [R, SumR, error, dt_A] = euler_solver(method,V,E,inedges,bedges,ufs,q,R,S,CFL,gamma);
    
    Rnorm(N) = sqrt(SumR);

    if  Rnorm(N) < 10^(-5)
        break
    end
    
    u_safe = q;

    % Check CFL for errors
    if error == 1 
        fprintf('Error occured!!: Back off on CFL! \n')
        CFLTuner = 1; % Back off on the CFL
        q = u_safe; % rReset to previous safe state
        continue 
    end  

    CFL = CFLmax - (0.5*CFLTuner);  % maybe this helps??

    q = qo - dt_A.*R; 

    % Possible 2nd Stage: (maybe help error problem??)
    %[R, SumR, error, dt_A] = euler_solver(method,V,E,inedges,bedges,ufs,q,R,S,e_ind,CFL,gamma); 

    %q = 0.75*qo + 0.25*(q - dt_A*R);

    dt = dt_A;
    
    if count >= count_max
        outputstate = fopen('U.txt','w');
        fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E\n',q');
        fclose(outputstate);
    
        count = 0;
        % Output the current counter: aka N = time passed and t being the
        % current
        fprintf('Working: (N = %d & t = %d & dt = %d & dt_A = %d)\n', N, t, dt, dt_A);
        
        CFLTuner = 0; % Reset tuner
    end  

    N = N + 1;
end

if N == Niter
    fprintf('/n...CAUTION...\n\n');
end

fprintf('/n...time is \n\n');
cputime = toc; disp(['CPU time: ' num2str(cputime),' s'])

%%

% Calculate converged flow properties
r = q(:,1);
U = 1./q(:,1).*sqrt(q(:,2).^2+q(:,3).^2); 
P = (gamma-1).*(q(:,4)-0.5.*q(:,1).*U.^2); 
a = (gamma*P./q(:,1)).^(0.5);
M = U./a; 
e = P./((gamma-1)*q(:,1));   
ss = log(P./r.^gamma); 

x = V(:,1);
y = V(:,2);

fprintf('<<<<< EULER SOLVER COMPLETE >>>>>\n\n');

save = input('Do you want to save the data? (y/n)', 's'); %maybe ask user for file name, if i have time
if (strcmp(save,'y') || strcmp(save,'yes') || strcmp(save,'Yes'))
    switch save
        case 'yes'
            outputstate = fopen('U_roe_40nx_discretize60.txt','w');
            fprintf(outputstate,'%23.16E %23.16E %23.16E %23.16E\n',q');
            fclose(outputstate);
    end
end

%% Plot outputs
fig = 0;
fig = fig + 1;
figure(fig)
clf(fig)
meshplot(M,E,V,'Local Mach Number','Converged Mach Field',1,0);
hold on 
set(findall(gcf,'-property','FontSize'),'FontSize',12);

%%

%
fig = fig + 1;
figure(fig)
clf(fig)
meshplot(P,E,V,'Local Pressure','Converged Pressure Field',1,0)
hold on 
set(findall(gcf,'-property','FontSize'),'FontSize',12);
%saveas(gcf,strcat(saveprefix,'ConvPress.jpg'));


%%
fig = fig + 1;
figure(fig)
clf(fig)
meshplot(M,E,V,'Local Mach Number','Converged Mach Field',1,0);
xlim([0, (x_ramp_start+x_length)])
ylim([0, 1])
hold on 
set(findall(gcf,'-property','FontSize'),'FontSize',12);
%saveas(gcf,strcat(saveprefix,'ConvMachZoom.jpg'));


%% Plot Mesh
answer = input('Would you like plot the mesh? (y/n): ','s');

if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes'))
    fig = fig + 1;
    figure(fig)
    clf(fig)
    hold on
    nulldat = zeros(size(E,1));
    meshplot(nulldat,E,V,'colorlabel','Full Mesh',1,1)
    colorbar('off');
    %saveas(gcf,strcat(saveprefix,'MeshFull.jpg'))
end  

%%
answer = input('Would you like to plot any properites on bottom surface? (y/n): ','s');
if (strcmp(answer,'y') || strcmp(answer,'yes') || strcmp(answer,'Yes'))
    
    prop_in = 'Pressure';
    prop = P;
    
    Mfs = 3;
    
    ax0 = [min(V(:,1)) max(V(:,1)) min(V(:,2)) max(V(:,2))]; 
    
    limx_neg = ax0(1)*1.75;
    limx_pos = ax0(2)*1.75;
    limy_neg = ax0(3)*1.75;
    limy_pos = ax0(4)*1.75;
    
    datamatrix = [];
    
    for i = (1:length(E(:,1)))
        
        tri = zeros(3,2); 
        
        p1 = E(i,1); 
        p2 = E(i,2);
        p3 = E(i,3); 
         
        tri(1,1) = V(p1,1); 
        tri(2,1) = V(p2,1);
        tri(3,1) = V(p3,1);
     
        tri(1,2) = V(p1,2); 
        tri(2,2) = V(p2,2);
        tri(3,2) = V(p3,2);
       
        if min(tri(:,1)) >= limx_neg && max(tri(:,1)) <= limx_pos && min(tri(:,2)) >= limy_neg && max(tri(:,2)) <= limy_pos
            datamatrix = [datamatrix; tri(:,1)', tri(:,2)', prop(i)]; % append values to matrix
        end
    end
    
    extract_line = [datamatrix(1,1) datamatrix(1,4) datamatrix(1,7)];
    j = 2;
    for i = 1:(n_y-1)
        x_values = datamatrix(j,1);
        y_values= datamatrix(j,4);
        values= datamatrix(j,7);
        extract_line = [extract_line; x_values, y_values, values];
        j = j + 2;
    end
    
    if thetad == 16.6992
        P2_P1 = 6.68634592;
    elseif thetad == 30.9638
        P2_P1 = 3.12303510;
    else
        P2_P1 = 2.19465313;
    end
    
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
