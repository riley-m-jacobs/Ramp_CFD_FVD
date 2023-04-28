%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Euler Inviscid Compressible Flow Solver
% Compression Ramp Mesh w/ Roe, HLLE Flux Formulas
%
% Riley Jacobs
% 04.26.2023
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close
clc

% Graph headers
out.header = 'M = 2, thetad = 15';
out.datasaver = 'Nx_20_Disc_10_M_2_theta_15';

% Characterize the problem you want:
out.Mfs = 2;              % Mach number of incoming free stream flow
out.thetad = 15;          % Turning angle of ramp in degrees
out.alphafs = 0;          % Enter theta direction of fs, typically 0 degrees
out.gamma = 1.4;          % Specific heat ratio

out.n_x = 20;             % Discretize the pre-ramp region, 3x in post-ramp region
out.discretize = 20;      % Discretize the number of nodes up the ramp
out.method = 'Roe';      % Roe or HLLE
out.CFLmax = 0.5;         % CFL condition --> max starting
out.sanity = false;

% Start from converged state?
out.state = false;

% Plotting 
out.data_save = false;             % Save data
out.mesh_plot = false;             % Plot the empty mesh
out.mach_plot = true;              % Plot the mach contours
out.pres_plot = true;              % Plot the pressure contours
out.save_plots = false;            % Save the contours
out.pres_bot = false;              % Plot the bottom pressure
out.mach_bot = true;              % Plot the mach along bottom
out.zoom_plot = false;             % Save the zoomed in wedge plot

% Run the CFD Ramp Solver
[sol, error] = driver(out);
