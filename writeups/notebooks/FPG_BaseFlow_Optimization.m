% Script file loading necessary packages for use within Octave
1; if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end

%% Define parameters for nozzle cases of interest
Ma = 1; gam0 = 1.4; Rinner = 1; Router = 2;

%% Solve subsonic nozzle (specifying inflow) and plot to file
nozzle(Ma, gam0, Router, Rinner, -2/7, 1, 1);  % Breaks for -1/3
title('Subsonic nozzle: outflow <- inflow');
print('nozzle_subsonic.eps', '-depsc2', '-S512,384', '-F:8');
close();

%% Subsonic nozzles may (more robustly) be defined with nearly sonic outflows
% nozzle(Ma, gam0, Rinner, Router, -1/Ma+sqrt(eps), 1, 1);

%% Solve supersonic nozzle (specifying inflow) and plot to file
nozzle(Ma, gam0, Rinner, Router,  1/Ma+sqrt(eps), 1, 1);
title('Supersonic nozzle: inflow -> outflow');
print('nozzle_supersonic.eps', '-depsc2', '-S512,384', '-F:8');
close();
