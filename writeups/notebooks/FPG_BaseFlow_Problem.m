% Script file loading necessary packages for use within Octave
1; if exist('OCTAVE_VERSION') ~= 0; pkg load odepkg; end

% Solve nozzle initial value problem for u(r) on [R1, R2] given u1(R1)
% via an ODE integrator-ready implementation of velocity derivative u'(r, u).
% Afterward, compute sound speed squared, density, and pressure from r, u.
function [r, u, rho, p, a2, up, pp] = nozzle(Ma, gam0, R1, R2, u1, rho1, p1)

  Ma2   = Ma*Ma;
  up    = @(r,u) -(u./r) .* (2 + Ma2*(gam0 - 1) - Ma2*(gam0 - 1)*u.**2) ...
                         ./ (2 + Ma2*(gam0 - 1) - Ma2*(gam0 + 1)*u.**2);
  tol   = sqrt(eps);
  vopt  = odeset('RelTol',      tol, 'AbsTol',  tol, ...
                 'InitialStep', tol, 'MaxStep', sqrt(sqrt(tol)));
  [r,u] = ode45(up, [R1 R2], u1, vopt);
  up    = up(r, u); % Shadow
  a2    = 1 + 0.5*Ma2*(gam0 - 1)*(1 - u.**2);
  rho   = rho1 * exp(-2*pi*Ma2 * cumtrapz(r, r.*u.*up ./ a2));
  pp    = - Ma2 * rho.*u.*up;
  p     = p1 + 2*pi * cumtrapz(r, r.*pp);

  if 0 == nargout
    figure();
    plot(r, u, 'r-', r, rho, 'b-', r, p, 'g-');
    legend('Velocity', 'Density', 'Pressure', 'location', 'northwest');
    xlabel('Radius');
  end

endfunction

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
